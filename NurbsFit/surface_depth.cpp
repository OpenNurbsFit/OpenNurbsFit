/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2012-, Thomas Mörwald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 */

#include "surface_depth.h"
#include <stdexcept>

//#include <glpk.h>
//#include <gurobi/gurobi_c++.h>

using namespace onf;

FittingSurfaceDepth::ROI::ROI(const Eigen::MatrixXd &points)
{
  double x_min(DBL_MAX), x_max(DBL_MIN);
  double y_min(DBL_MAX), y_max(DBL_MIN);

  for(Eigen::MatrixXd::Index i=0; i<points.rows(); i++)
  {
    const Eigen::Vector3d& p = points.row(i);

    if(p(0)<x_min)
      x_min=p(0);
    if(p(0)>x_max)
      x_max=p(0);

    if(p(1)<y_min)
      y_min=p(1);
    if(p(1)>y_max)
      y_max=p(1);
  }

  x = x_min;
  y = y_min;
  width = x_max-x_min;
  height = y_max-y_min;
}

FittingSurfaceDepth::FittingSurfaceDepth(int order,
                                         int cps_width, int cps_height,
                                         ROI img_roi,
                                         const Eigen::MatrixXd& points)
  : m_quiet(true), m_solver(NULL)
{
  if(points.cols()!=3)
    throw std::runtime_error("[FittingSurfaceDepth::FittingSurfaceDepth] Error, points must be a matrix Nx3.");

  initSurface(order, cps_width, cps_height, img_roi);
  initSolver(points);
  solve(points.col(2));
}

FittingSurfaceDepth::FittingSurfaceDepth(int order,
                                         int cps_width, int cps_height,
                                         ROI img_roi,
                                         const Eigen::MatrixXd& points,
                                         const std::vector<int>& indices)
  : m_quiet(true), m_solver(NULL)
{
  if(points.cols()!=3)
    throw std::runtime_error("[FittingSurfaceDepth::FittingSurfaceDepth] Error, points must be a matrix Nx3.");

  initSurface(order, cps_width, cps_height, img_roi);
  initSolver(points, indices);
  solve(points.col(2), indices);
}

FittingSurfaceDepth::~FittingSurfaceDepth()
{
  if(m_solver!=NULL)
    delete m_solver;
}


void FittingSurfaceDepth::initSurface(int order, int cps_width, int cps_height, ROI img_roi)
{
  if(cps_width<order)
    cps_width=order;
  if(cps_height<order)
    cps_height=order;

  m_roi = img_roi;
  m_nurbs = ON_NurbsSurface (1, false, order, order, cps_width, cps_height);

  double ddx = m_roi.width  / (m_nurbs.KnotCount(0) - 2*(order-2) - 1);
  double ddy = m_roi.height / (m_nurbs.KnotCount(1) - 2*(order-2) - 1);

  m_nurbs.MakeClampedUniformKnotVector (0, ddx);
  m_nurbs.MakeClampedUniformKnotVector (1, ddy);

  for (int i = 0; i < m_nurbs.KnotCount(0); i++)
  {
    double k = m_nurbs.Knot (0, i);
    m_nurbs.SetKnot (0, i, k + m_roi.x);
  }

  for (int i = 0; i < m_nurbs.KnotCount(1); i++)
  {
    double k = m_nurbs.Knot (1, i);
    m_nurbs.SetKnot (1, i, k + m_roi.y);
  }

  m_b = Eigen::VectorXd(m_nurbs.CVCount(),1);
  for (int i = 0; i < m_nurbs.CVCount(0); i++)
  {
    for (int j = 0; j < m_nurbs.CVCount(1); j++)
    {
      m_nurbs.SetCV (i, j, ON_3dPoint(0,0,0));
      m_b(grc2gl(i,j),0) = 0.0;
    }
  }

//    ON_TextLog out;
//    m_nurbs.Dump (out);
}

void FittingSurfaceDepth::initSolver(const Eigen::MatrixXd &points)
{
  if(points.cols()!=3)
    throw std::runtime_error("[FittingSurfaceDepth::initSolver] Error, points must be a matrix Nx3.");

  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FittingSurfaceDepth::initSolver] Error, surface not initialized (initSurface).");

  m_K = SparseMatrix( points.rows(), m_nurbs.CVCount() );

  typedef Eigen::Triplet<double> Tri;
  std::vector<Tri> tripletList;
  tripletList.resize( points.rows() * m_nurbs.Order(0) * m_nurbs.Order(1) );

  if(!m_quiet)
    printf("[FittingSurfaceDepth::initSolver] entries: %lu  rows: %lu  cols: %d\n",
           points.rows()* m_nurbs.Order(0) * m_nurbs.Order(1),
           points.rows(),
           m_nurbs.CVCount());

  double *N0 = new double[m_nurbs.Order (0) * m_nurbs.Order (0)];
  double *N1 = new double[m_nurbs.Order (1) * m_nurbs.Order (1)];
  int E,F;
  int ti(0);

  for(Eigen::MatrixXd::Index row=0; row<points.rows(); row++)
  {
    const Eigen::Vector3d& p = points.row(row);

    E = ON_NurbsSpanIndex (m_nurbs.m_order[0], m_nurbs.m_cv_count[0], m_nurbs.m_knot[0], p(0), 0, 0);
    F = ON_NurbsSpanIndex (m_nurbs.m_order[1], m_nurbs.m_cv_count[1], m_nurbs.m_knot[1], p(1), 0, 0);

    ON_EvaluateNurbsBasis (m_nurbs.Order (0), m_nurbs.m_knot[0] + E, p(0), N0);
    ON_EvaluateNurbsBasis (m_nurbs.Order (1), m_nurbs.m_knot[1] + F, p(1), N1);

    for (int i = 0; i < m_nurbs.Order (0); i++)
    {
      for (int j = 0; j < m_nurbs.Order (1); j++)
      {
        tripletList[ti] = Tri( row, lrc2gl (E, F, i, j), N0[i] * N1[j] );
        ti++;
      } // j
    } // i

  } // row

  delete [] N1;
  delete [] N0;

  m_K.setFromTriplets(tripletList.begin(), tripletList.end());

  if(m_solver==NULL)
    m_solver = new SPQR();
  m_solver->compute(m_K);
  if(m_solver->info()!=Eigen::Success)
    throw std::runtime_error("[FittingSurfaceDepth::initSolver] decomposition failed.");

  if(!m_quiet)
    printf("[FittingSurfaceDepth::initSolver] decomposition done\n");

  m_use_indices = false;
}

void FittingSurfaceDepth::initSolver(const Eigen::MatrixXd &points, const std::vector<int>& indices)
{
  if(points.cols()!=3)
    throw std::runtime_error("[FittingSurfaceDepth::initSolver] Error, points must be a matrix Nx3.");
  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FittingSurfaceDepth::initSolver] Error, surface not initialized (initSurface).");

  //  m_nsolver.assign(m_roi.width*m_roi.height, m_nurbs.CVCount(), 1);

  m_K = SparseMatrix( indices.size(), m_nurbs.CVCount() );

  typedef Eigen::Triplet<double> Tri;
  std::vector<Tri> tripletList;
  tripletList.resize( indices.size() * m_nurbs.Order(0) * m_nurbs.Order(1) );

  if(!m_quiet)
    printf("[FittingSurfaceDepth::initSolver] entries: %lu  rows: %lu  cols: %d\n",
           indices.size() * m_nurbs.Order(0) * m_nurbs.Order(1),
           indices.size(),
           m_nurbs.CVCount());

  double *N0 = new double[m_nurbs.Order (0) * m_nurbs.Order (0)];
  double *N1 = new double[m_nurbs.Order (1) * m_nurbs.Order (1)];
  int E,F;
  int ti(0);

  for(size_t row=0; row<indices.size(); row++)
  {
    const Eigen::Vector3d& p = points.row(indices[row]);

    E = ON_NurbsSpanIndex (m_nurbs.m_order[0], m_nurbs.m_cv_count[0], m_nurbs.m_knot[0], p(0), 0, 0);
    F = ON_NurbsSpanIndex (m_nurbs.m_order[1], m_nurbs.m_cv_count[1], m_nurbs.m_knot[1], p(1), 0, 0);

    ON_EvaluateNurbsBasis (m_nurbs.Order (0), m_nurbs.m_knot[0] + E, p(0), N0);
    ON_EvaluateNurbsBasis (m_nurbs.Order (1), m_nurbs.m_knot[1] + F, p(1), N1);

    for (int i = 0; i < m_nurbs.Order (0); i++)
    {
      for (int j = 0; j < m_nurbs.Order (1); j++)
      {
        tripletList[ti] = Tri( row, lrc2gl (E, F, i, j), N0[i] * N1[j] );
        //          m_nsolver.K(ti,lrc2gl (E, F, i, j),N0[i]*N1[j]);
        ti++;
      } // j
    } // i

  } // row

  delete [] N1;
  delete [] N0;

  m_K.setFromTriplets(tripletList.begin(), tripletList.end());

  if(m_solver==NULL)
    m_solver = new SPQR();
  m_solver->compute(m_K);
  if(m_solver->info()!=Eigen::Success)
    throw std::runtime_error("[FittingSurfaceDepth::initSolver] decomposition failed.");
  if(!m_quiet)
    printf("[FittingSurfaceDepth::initSolver] decomposition done\n");
  m_use_indices = true;
}

void FittingSurfaceDepth::solve(const Eigen::VectorXd &z)
{
  if(m_use_indices)
    throw std::runtime_error("[FittingSurfaceDepth::solve] Error, solver initialized with indices (use solve(Eigen::VectorXd&, const std::vector<int>&) instead.\n");

  m_b = m_solver->solve(z);

  updateSurf();
}

void FittingSurfaceDepth::solve(const Eigen::VectorXd &z, const std::vector<int>& indices)
{
  if(!m_use_indices)
    throw std::runtime_error("[FittingSurfaceDepth::solve] Error, solver initialized without indices (use solve(Eigen::VectorXd&) instead.\n");

  Eigen::VectorXd zn(m_solver->rows(),1);

  for(size_t i=0; i<indices.size(); i++)
    zn(i) = z(indices[i]);

  m_b = m_solver->solve(zn);

  updateSurf();
}

void FittingSurfaceDepth::updateSurf()
{
  int ncp = m_nurbs.CVCount ();

  if(m_b.rows()!=ncp)
    throw std::runtime_error("[FittingSurfaceDepth::updateSurf] Error, number of control points does not match.");

  for (int i = 0; i < ncp; i++)
  {
    ON_3dPoint cp;

    cp.x = m_b(i, 0);
    cp.y = cp.x;
    cp.z = cp.x;

    m_nurbs.SetCV (gl2gr(i), gl2gc(i), cp);
  }
}

Eigen::VectorXd FittingSurfaceDepth::GetError(const Eigen::VectorXd& z)
{
  if(m_use_indices)
    throw std::runtime_error("[FittingSurfaceDepth::GetError] Error, solver initialized with indices (use solve(Eigen::VectorXd&, const std::vector<int>&) instead.\n");

  // compute K*b (i.e. points on surface)
  Eigen::VectorXd s(z.rows(),1);
  s.setZero();
  for (int k=0; k<m_K.outerSize(); ++k)
  for (SparseMatrix::InnerIterator it(m_K,k); it; ++it)
    s(it.row()) += it.value() * m_b(it.col());

  // return (K*b-z)
  return (s-z);
}

Eigen::VectorXd FittingSurfaceDepth::GetError(const Eigen::VectorXd& z, const std::vector<int>& indices)
{
  if(!m_use_indices)
    throw std::runtime_error("[FittingSurfaceDepth::GetError] Error, solver initialized without indices (use solve(Eigen::VectorXd&) instead.\n");

  // compute K*b (i.e. points on surface)
  Eigen::VectorXd e(indices.size(),1);
  e.setZero();
  for (int k=0; k<m_K.outerSize(); ++k)
  for (SparseMatrix::InnerIterator it(m_K,k); it; ++it)
    e(it.row()) += it.value() * m_b(it.col());

  // compute (K*b-z)
  for(size_t i=0; i<indices.size(); i++)
    e(i) -= z(indices[i]);

  return e;
}

void onf::IncreaseDimension( const ON_NurbsSurface& src, ON_NurbsSurface& dest, int dim )
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
