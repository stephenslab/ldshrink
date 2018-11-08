#pragma once

#include <RcppEigen.h>
#include <limits>
#include <memory>
#include <mutex>

#include<thread>

class LDshrinkCor{
  const double m;
  const double ne;
  const double cutoff;
  const double GenoMult;
  const double theta;
  const double pre_mult;
  using	datapair=std::shared_ptr<std::pair<const Eigen::VectorXd,double> >;
public:
  LDshrinkCor(const double m_, const double ne_,const double cutoff_,const bool isGenotype=true):
    m(m_),
    ne(ne_),
    cutoff(cutoff_),
    GenoMult(isGenotype ?	0.5 : 1),
    theta(LDshrinkCor::calc_theta(m))
    ,pre_mult(        GenoMult * (1 - theta) * (1 - theta) + 0.5 * theta * (1 - 0.5 * theta))
  {
  }
  double check(double map_dist)const{
    double rho = 4 * ne * (map_dist) / 100;
    rho = -rho / (2 * m);
    rho = std::exp(rho);
    return( rho >= cutoff ? rho : 0);
  }
  double cor(const double rho,std::pair<datapair,datapair> inputs)const{
    const auto& x_1 =inputs.first->first;
    const double Nm1 = x_1.size()-1;
    const double& var_1 =inputs.first->second;
    const auto& x_2 =inputs.second->first;
    const double& var_2 =inputs.second->second;
    double r= ((x_1.dot(x_2)/Nm1)*rho*GenoMult*(1-theta)*(1-theta));
    r=r/(pre_mult*std::sqrt(var_1*var_2));
    return(r);
  };
private:
  static double calc_nmsum(const double m) {
  int msize = (2 * (int)m - 1);
  Eigen::ArrayXd tx(msize);
  tx.setLinSpaced(msize, 1, (int)(2 * m - 1));
  return (1 / tx).sum();
  }

  static double calc_theta(const double m){
    double nmsum=calc_nmsum(m);
    return((1/nmsum)/(2*m+1/nmsum));
  }


};


class SampleCor{
  const double anno_cutoff;
  const double cor_cutoff;
  using	datapair=std::shared_ptr<std::pair<const Eigen::VectorXd,double> >;


public:
  SampleCor(const double anno_cutoff_ = std::numeric_limits<double>::max(),
            const double cor_cutoff_ = 0)
    : anno_cutoff(anno_cutoff_), cor_cutoff(cor_cutoff_*cor_cutoff) //,
                                                           // m()
  {
    // Rcpp::Rcerr < "will kill all distances that are less than : " << anno_cutoff
    //                                                               << std::endl;
  }
  double check(const double anno_dist)const{
   // std::lock_guard<std::mutex> myLock(m);
  //  Rcpp::Rcerr<<anno_dist<<std::endl;
    return(anno_dist<=anno_cutoff ? anno_dist : 0);
  }
  double cor(const double rho,std::pair<datapair,datapair> inputs)const{
    const auto& x_1 =inputs.first->first;
    const double Nm1 = x_1.size()-1;
    const double& var_1 =inputs.first->second;
    const auto& x_2 =inputs.second->first;
    const double& var_2 =inputs.second->second;
    const double tcor =	((x_1.array()/std::sqrt(var_1)).matrix().dot((x_2.array()/std::sqrt(var_2)).matrix()))/Nm1;
    return tcor*tcor >= cor_cutoff ? tcor : 0;
  };
};
