"""Exercise the real shared compass renderer with valid and unavailable headings."""
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class PiiDeviceHeadingTest(unittest.TestCase):
    def test_invalid_compass_does_not_report_north(self):
        stub = r'''
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
constexpr int TFT_BLACK=0,TFT_DARKGREY=1,TFT_LIGHTGREY=2,TFT_WHITE=3,
              TFT_RED=4,TFT_GREEN=5,TFT_ORANGE=6;
inline std::vector<std::string> printed;
inline int needles=0;
struct Screen {int width(){return 128;}int height(){return 128;}};
inline struct {Screen Display;} M5;
struct M5Canvas {
  std::vector<uint8_t> buffer;
  void setColorDepth(int){}
  void* createSprite(int w,int h){buffer.resize(unsigned(w*h*2));return buffer.data();}
  void deleteSprite(){buffer.clear();}
  void* getBuffer(){return buffer.data();}
  void pushSprite(Screen*,int,int){}
  void fillSprite(int){}
  void drawCircle(int,int,int,int){}
  void drawLine(int,int,int,int,int){}
  void fillTriangle(int,int,int,int,int,int,int){++needles;}
  void fillCircle(int,int,int,int){}
  void fillRect(int,int,int,int,int){}
  void setTextColor(int,int){}
  void setTextSize(int){}
  int textWidth(const char* s){return int(std::strlen(s))*6;}
  int fontHeight(){return 8;}
  void setCursor(int,int){}
  void print(const char* s){printed.emplace_back(s);}
};
'''
        cpp = r'''
#include "AtomS3R/AtomS3R_CompassUI.h"
#include <algorithm>
#include <cassert>
#include <limits>
bool text(const char* s){return std::find(printed.begin(),printed.end(),s)!=printed.end();}
int main(){
  CompassUI ui;ui.begin();assert(ui.ok());
  printed.clear();needles=0;ui.draw(137,true,51,false);
  assert(needles==1&&text("137")&&!text("WAIT"));
  for(float heading:{0.f,137.f,std::numeric_limits<float>::quiet_NaN(),
                     std::numeric_limits<float>::infinity()}){
    printed.clear();needles=0;ui.draw(heading,false,51,false);
    assert(needles==0&&text("---")&&text("WAIT")&&!text("0")&&!text("N"));
  }
  printed.clear();needles=0;ui.draw(std::numeric_limits<float>::quiet_NaN(),true);
  assert(needles==0&&text("WAIT"));
  printed.clear();needles=0;ui.draw(0,true);
  assert(needles==1&&text("0")&&text("N"));
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / 'M5Unified.h').write_text(stub)
            source = folder / 'compass.cpp'
            binary = folder / 'compass'
            source.write_text(cpp)
            subprocess.run(['g++', '-std=c++20', '-Wall', '-Wextra', '-Werror',
                            '-I'+str(folder), '-I'+str(ROOT/'src'), str(source),
                            '-o', str(binary)], check=True, capture_output=True, text=True)
            subprocess.run([str(binary)], check=True, timeout=5)
