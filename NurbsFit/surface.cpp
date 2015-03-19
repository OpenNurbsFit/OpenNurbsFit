
#include "surface.h"
#include <stdexcept>

using namespace nurbsfit;

FitSurface::Domain::Domain(const Eigen::VectorXd& param0, const Eigen::VectorXd& param1)
{
  if(param0.rows()!=param1.rows())
    throw std::runtime_error("[FitSurface::Domain::Domain] Error, param vectors must be of same length.");

  double x_min(DBL_MAX), x_max(DBL_MIN);
  double y_min(DBL_MAX), y_max(DBL_MIN);

  for(Eigen::MatrixXd::Index i=0; i<param0.rows(); i++)
  {
    const double& p0 = param0(i);
    const double& p1 = param1(i);

    if(p0<x_min)
      x_min=p0;
    if(p0>x_max)
      x_max=p0;

    if(p1<y_min)
      y_min=p1;
    if(p1>y_max)
      y_max=p1;
  }

  x = x_min;
  y = y_min;
  width = x_max-x_min;
  height = y_max-y_min;
  border_offset = 0;
}

FitSurface::~FitSurface()
{
  if(m_solver!=NULL)
    delete m_solver;
}

void nurbsfit::IncreaseDimension( const ON_NurbsSurface& src, ON_NurbsSurface& dest, int dim )
{
  dest.m_dim          = dim;
  dest.m_is_rat       = src.m_is_rat;
  dest.m_order[0]     = src.m_order[0];
  dest.m_order[1]     = src.m_order[1];
  dest.m_cv_count[0]  = src.m_cv_count[0];
  dest.m_cv_count[1]  = src.m_cv_count[1];
  dest.m_cv_stride[1] = dest.m_is_rat ? dest.m_dim+1 : dest.m_dim;
  dest.m_cv_stride[0] = dest.m_cv_count[1]*dest.m_cv_stride[1];
  if ( src.m_knot[0] )
  {
    // copy knot array
    dest.ReserveKnotCapacity( 0, dest.KnotCount(0) );
    memcpy( dest.m_knot[0], src.m_knot[0], dest.KnotCount(0)*sizeof(*dest.m_knot[0]) );
  }
  if ( src.m_knot[1] )
  {
    // copy knot array
    dest.ReserveKnotCapacity( 1, dest.KnotCount(1) );
    memcpy( dest.m_knot[1], src.m_knot[1], dest.KnotCount(1)*sizeof(*dest.m_knot[1]) );
  }
  if ( src.m_cv )
  {
    // copy cv array
    dest.ReserveCVCapacity( dest.m_cv_count[0]*dest.m_cv_count[1]*dest.m_cv_stride[1] );
    const int dst_cv_size = dest.CVSize()*sizeof(*dest.m_cv);
    const int src_stride[2] = {src.m_cv_stride[0],src.m_cv_stride[1]};
    if ( src_stride[0] == dest.m_cv_stride[0] && src_stride[1] == dest.m_cv_stride[1] )
    {
      memcpy( dest.m_cv, src.m_cv, dest.m_cv_count[0]*dest.m_cv_count[1]*dest.m_cv_stride[1]*sizeof(*dest.m_cv) );
    }
    else
    {
      const double *src_cv;
      double *dst_cv = dest.m_cv;
      int i, j;
      for ( i = 0; i < dest.m_cv_count[0]; i++ )
      {
        src_cv = src.CV(i,0);
        for ( j = 0; j < dest.m_cv_count[1]; j++ )
        {
          memcpy( dst_cv, src_cv, dst_cv_size );
          dst_cv += dest.m_cv_stride[1];
          src_cv += src_stride[1];
        }
      }
    }
  }
}
