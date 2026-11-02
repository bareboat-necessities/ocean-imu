#!/usr/bin/env python3
"""Exercise the real PSRAM preflight before any full-calibration captures."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "src/AtomS3R/AtomS3R_ImuCalWizard.h").read_text()


def block(marker):
    start = source.index(marker)
    brace = source.index("{", start)
    depth, end = 1, brace + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[start:end]


flow = block("bool runAndSave_(")
assert "prepareGyroRotationStorage_()" in flow, "PSRAM must be reserved before full captures"
assert (flow.index("StartAction::NOISE") < flow.index("prepareGyroRotationStorage_()")
        < flow.index("accel_.runMainStage")), "noise must stay separate from full-capture PSRAM admission"
capture = block("bool captureGyroRotations_(")
assert "gyro_rotations_.reset()" not in capture, "capture must retain the preflight reservation"

cpp = r'''
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <new>
#include <string>
#include <vector>
#include <iostream>
constexpr unsigned MALLOC_CAP_SPIRAM=1,MALLOC_CAP_8BIT=2;
size_t total=0,available=0,largest=0,allocated_bytes=0;
unsigned allocations=0,frees=0,last_caps=0;
bool found=false,fail_allocate=false;
bool psramFound(){return found;}
size_t heap_caps_get_total_size(unsigned){return total;}
size_t heap_caps_get_free_size(unsigned){return available;}
size_t heap_caps_get_largest_free_block(unsigned){return largest;}
void* heap_caps_malloc(size_t bytes,unsigned caps) {
  ++allocations;last_caps=caps;allocated_bytes=bytes;
  return fail_allocate?nullptr:std::malloc(bytes);
}
void heap_caps_free(void* p){++frees;std::free(p);}
struct {
  std::string messages;
  void println(const char* s){messages+=s;messages+='\n';}
  template<class... A> void printf(const char* fmt,A... args){
    char text[512];std::snprintf(text,sizeof(text),fmt,args...);messages+=text;
  }
} Serial;
struct RotationFit {char samples[1024]{};};
void require(bool b,const char* m){if(!b){std::cerr<<m<<'\n';std::exit(1);}}
class ImuCalWizard {
public:
  struct UI {
    std::vector<std::string> messages;
    void failLines(const char* title,const char* a,const char* b=nullptr,const char* c=nullptr){
      for(const char* s:{title,a,b,c})if(s)messages.emplace_back(s);
    }
    void fail(const char* title,const char* reason){failLines(title,reason);}
  } ui_;
DELETER;
  std::unique_ptr<RotationFit,RotationDeleter> gyro_rotations_;
PREPARE
};
int main(){
  {
    ImuCalWizard w;
    require(!w.prepareGyroRotationStorage_() && allocations==0 && !w.gyro_rotations_,
            "unavailable PSRAM must fail before allocating or capturing");
    require(w.ui_.messages.size()>=3 && Serial.messages.find("PSRAM")!=std::string::npos,
            "missing PSRAM must explain board configuration");
    bool opi=false;for(const auto& s:w.ui_.messages){
      opi=opi || s.find("OPI")!=std::string::npos;
      require(s.size()<=20,"preflight UI text exceeds display width");
    }
    require(opi,"missing PSRAM must request OPI mode");
  }
  total=8*1024*1024;available=total;largest=sizeof(RotationFit)-1;found=true;
  {
    ImuCalWizard w;
    require(!w.prepareGyroRotationStorage_() && allocations==0,
            "fragmentation must be distinguished from total free bytes");
  }
  largest=total;fail_allocate=true;
  {
    ImuCalWizard w;
    require(!w.prepareGyroRotationStorage_() && allocations==1 && !w.gyro_rotations_,
            "allocation failure after preflight must fail closed");
  }
  fail_allocate=false;found=false; // registered SDK heap is authoritative
  {
    ImuCalWizard w;
    require(w.prepareGyroRotationStorage_() && w.gyro_rotations_ && allocations==2,
            "registered external heap must support the reservation");
    auto* reserved=w.gyro_rotations_.get();
    require(last_caps==(MALLOC_CAP_SPIRAM|MALLOC_CAP_8BIT) && allocated_bytes==sizeof(RotationFit),
            "rotation trajectory must stay in byte-addressable external RAM");
    require(w.prepareGyroRotationStorage_() && w.gyro_rotations_.get()==reserved && allocations==2,
            "capture must reuse the reservation without a second allocation");
  }
  require(frees==1,"reservation must be freed exactly once");
  std::cout<<"wizard PSRAM preflight: unavailable/fragmented/allocation failure/reservation/reuse PASS\n";
}
'''
cpp = (cpp.replace("DELETER", block("struct RotationDeleter"))
       .replace("PREPARE", block("bool prepareGyroRotationStorage_(")))
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / "wizard-memory"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-x", "c++", "-", "-o", str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True, timeout=10)
