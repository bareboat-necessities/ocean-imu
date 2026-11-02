#define EIGEN_NON_ARDUINO
#include "imu_calibrate/SensorNoise.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
// Replay [NOISERAW] lines using exactly the MCU acquisition and identification
// code. No second fitting implementation or new Python dependency.
int main(int argc,char** argv) {
  if(argc!=2){std::cerr<<"usage: sensor_noise-replay capture.log\n";return 2;}
  std::ifstream in(argv[1]);if(!in){std::cerr<<"cannot read capture\n";return 2;}
  imu_cal::StationaryNoiseCapture<> capture;bool begun=false;std::string line;unsigned lines=0;
  while(std::getline(in,line)) {
    const auto prefix=line.find("[NOISERAW] ");if(prefix==std::string::npos)continue;
    std::string data=line.substr(prefix+11);std::replace(data.begin(),data.end(),',',' ');
    std::istringstream row(data);uint32_t us;float ax,ay,az,wx,wy,wz,temp;std::string extra;
    if(!(row>>us>>ax>>ay>>az>>wx>>wy>>wz>>temp) || row>>extra){std::cerr<<"malformed raw record\n";return 1;}
    if(!begun){capture.begin(us);begun=true;}
    const auto status=capture.add(us,Eigen::Vector3f(ax,ay,az),Eigen::Vector3f(wx,wy,wz),temp);
    ++lines;
    if(status!=imu_cal::NoiseCaptureStatus::CAPTURING && status!=imu_cal::NoiseCaptureStatus::READY) {
      std::cerr<<"capture rejected at row "<<lines<<", status "<<unsigned(status)<<"\n";return 1;
    }
  }
  imu_cal::SensorNoiseModel model;
  if(!capture.report(model)){std::cerr<<"40-minute capture or identifiable white region missing\n";return 1;}
  std::cout<<std::setprecision(9)<<"{\"duration_ms\":"<<model.duration_ms<<",\"samples\":"<<model.samples
    <<",\"period_s\":"<<model.period_s<<",\"white_density\":[";
  for(int i=0;i<6;++i)std::cout<<(i?",":"")<<model.white_density[i];
  std::cout<<"],\"bias_mask\":"<<unsigned(model.bias_mask)<<",\"walk_mask\":"<<unsigned(model.walk_mask)
    <<",\"allan\":[";
  for(int l=0;l<17;++l) {
    std::cout<<(l?",":"")<<"{\"tau\":"<<capture.tau(l)<<",\"pairs\":"<<capture.pairs(l)<<",\"adev\":[";
    for(int i=0;i<6;++i)std::cout<<(i?",":"")<<capture.deviation(l,i);
    std::cout<<"]}";
  }
  std::cout<<"]}\n";return 0;
}
