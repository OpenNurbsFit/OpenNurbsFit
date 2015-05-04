
#include "curve.h"
#include <stdexcept>

using namespace nurbsfit;

FitCurve::Domain::Domain(const Eigen::VectorXd &param)
{
  if(param.rows()<2)
    throw std::runtime_error("[FitCurve::Domain::Domain] Error, insufficient parameter points (<2).");

  double param_min(DBL_MAX);
  double param_max(DBL_MIN);

  for(Eigen::MatrixXd::Index i=0; i<param.rows(); i++)
  {
    const double& p = param(i);

    if(p<param_min)
      param_min=p;
    if(p>param_max)
      param_max=p;
  }

  start = param_min;
  end = param_max;
}
