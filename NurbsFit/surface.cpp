
#include "surface.h"

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
