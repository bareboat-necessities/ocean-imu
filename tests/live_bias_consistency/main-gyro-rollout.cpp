#ifndef EIGEN_NON_ARDUINO
#define EIGEN_NON_ARDUINO
#endif
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>
namespace gyro_rollout {int ou2(int,double,const std::string&);int ou3(int,double,const std::string&);int tfg(int,double,const std::string&);}
int main(int argc,char**argv)try{if(argc!=5)throw std::runtime_error("family seed duration output-prefix");std::string f=argv[1];int s=std::stoi(argv[2]);double t=std::stod(argv[3]);if(t<180||!std::isfinite(t))throw std::runtime_error("duration");return f=="OU2"?gyro_rollout::ou2(s,t,argv[4]):f=="OU3"?gyro_rollout::ou3(s,t,argv[4]):f=="TFG"?gyro_rollout::tfg(s,t,argv[4]):throw std::runtime_error("family");}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
