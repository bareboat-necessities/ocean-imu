#!/usr/bin/env python3
"""Exercise the actual wizard task dispatcher with a host scheduler stub."""
from pathlib import Path
import os
import re
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


stack = re.search(r"WIZARD_STACK_BYTES\s*=\s*(\d+);", source).group(1)
cpp = r'''
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
using TaskHandle_t=void*;using BaseType_t=int;
constexpr int pdPASS=1,ESP_OK=0;
struct ImuCalBlobV6 {unsigned token=0;};
struct ImuCalWizardCfg {static constexpr uint32_t WIZARD_STACK_BYTES=STACK;};
struct {void println(const char*){} template<class... A>void printf(const char*,A...){} } Serial;
std::thread worker;const auto caller=std::this_thread::get_id();
std::atomic<unsigned> feeds{0},deletes{0};
bool allocation_failure=false,watchdog_registered=true;unsigned requested_stack=0;int requested_core=-1;
int xPortGetCoreID(){return 1;}
void delay(unsigned){std::this_thread::yield();}
int esp_task_wdt_status(void*){return watchdog_registered?ESP_OK:-1;}
void esp_task_wdt_reset(){++feeds;}
void vTaskDelete(void*){++deletes;}
int xTaskCreatePinnedToCore(void(*fn)(void*),const char*,uint32_t bytes,void* ctx,
                           unsigned,TaskHandle_t* task,int core) {
  requested_stack=bytes;requested_core=core;
  if(allocation_failure){*task=nullptr;return 0;}
  *task=reinterpret_cast<void*>(1);worker=std::thread([=]{fn(ctx);});return pdPASS;
}
void require(bool value,const char* message){if(!value){std::cerr<<message<<'\n';std::exit(1);}}
class ImuCalWizard {
public:
  struct UI {unsigned failures=0;void fail(const char*,const char*){++failures;}} ui_;
  bool succeed=true,ran_off_caller=false;unsigned runs=0;
CONTEXT;
DISPATCH
CALLBACK
  static uint32_t hwmBytes_(){return requested_stack-8000;}
  bool runAndRestore_(ImuCalBlobV6& out) {
    ran_off_caller=std::this_thread::get_id()!=caller;
    if(watchdog_registered)while(!feeds.load())std::this_thread::yield();
    ++runs;if(succeed)out.token=42;return succeed;
  }
};
int main(){
  ImuCalWizard wizard;ImuCalBlobV6 out;
  require(wizard.runAndSave(out),"saved result was lost");worker.join();
  require(wizard.ran_off_caller && wizard.runs==1,"wizard body ran on the caller stack");
  require(requested_stack==ImuCalWizardCfg::WIZARD_STACK_BYTES && requested_stack>=32768,
          "wizard task stack must be passed in ESP-IDF bytes");
  require(requested_core==1 && feeds && deletes==1 && out.token==42,
          "core, watchdog, completion publication or output lifetime is wrong");
  wizard.succeed=false;out.token=7;feeds=0;
  require(!wizard.runAndSave(out),"cancelled wizard reported success");worker.join();
  require(out.token==7 && deletes==2,"cancellation changed the saved output");
  allocation_failure=true;
  require(!wizard.runAndSave(out) && wizard.ui_.failures==1 && wizard.runs==2 && out.token==7,
          "task allocation failure did not fail closed");
  allocation_failure=false;watchdog_registered=false;feeds=0;wizard.succeed=true;
  require(wizard.runAndSave(out),"unregistered caller watchdog blocked the wizard");worker.join();
  require(feeds==0 && deletes==3,"unregistered watchdog was reset");
  std::cout<<"wizard dispatcher: dedicated stack, caller core, save/cancel/allocation failure and watchdog PASS\n";
}
'''
cpp = (cpp.replace("=STACK;", "=" + stack + ";")
       .replace("CONTEXT", block("struct WizardTaskCtx"))
       .replace("DISPATCH", block("bool runAndSave("))
       .replace("CALLBACK", block("static void wizardTask_(")))
with tempfile.TemporaryDirectory() as directory:
    binary = Path(directory) / "wizard-task"
    subprocess.run(shlex.split(os.environ.get("CXX", "g++")) + [
        "-std=c++20", "-O1", "-pthread", "-x", "c++", "-", "-o", str(binary),
    ], input=cpp, text=True, check=True)
    subprocess.run([str(binary)], check=True, timeout=10)
