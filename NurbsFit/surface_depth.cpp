/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2015, Thomas Mörwald
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

using namespace nurbsfit;

FitSurfaceDepth::FitSurfaceDepth(int order0, int order1,
                                 int cps0, int cps1,
                                 Domain roi,
                                 const Eigen::MatrixXd& points)
{
  if(points.cols()!=3)
    throw std::runtime_error("[FitSurfaceDepth::FitSurfaceDepth] Error, points must be a matrix Nx3.");

  initSurface(order0, order1, cps0, cps1, roi);
  initSolver(points.col(0), points.col(1));
  solve(points.col(2));
}

FitSurfaceDepth::FitSurfaceDepth(int order0, int order1,
                                 int cps0, int cps1,
                                 Domain roi,
                                 const Eigen::MatrixXd& points,
                                 const std::vector<int>& indices)
{
  if(points.cols()!=3)
    throw std::runtime_error("[FitSurfaceDepth::FitSurfaceDepth] Error, points must be a matrix Nx3.");

  initSurface(order0, order1, cps0, cps1, roi);
  initSolver(points.col(0), points.col(1), indices);
  solve(points.col(2), indices);
}

void FitSurfaceDepth::initSurface(int order0, int order1, int cps0, int cps1, Domain roi)
{
  if( cps0 < order0 )
    cps0 = order0;

  if( cps1 < order1 )
    cps1 = order1;

  m_nurbs = ON_NurbsSurface (1, false, order0, order1, cps0, cps1);

  double ddx = roi.width  / (m_nurbs.KnotCount(0) - 2*(order0-2) - 1);
  double ddy = roi.height / (m_nurbs.KnotCount(1) - 2*(order1-2) - 1);

  m_nurbs.MakeClampedUniformKnotVector (0, ddx);
  m_nurbs.MakeClampedUniformKnotVector (1, ddy);

  for (int i = 0; i < m_nurbs.KnotCount(0); i++)
  {
    double k = m_nurbs.Knot (0, i);
    m_nurbs.SetKnot (0, i, k + roi.x);
  }

  for (int i = 0; i < m_nurbs.KnotCount(1); i++)
  {
    double k = m_nurbs.Knot (1, i);
    m_nurbs.SetKnot (1, i, k + roi.y);
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
}

void FitSurfaceDepth::initSolver(const Eigen::VectorXd& param0,
                                 const Eigen::VectorXd& param1)
{
  if(param0.rows()!=param1.rows())
    throw std::runtime_error("[FitSurfaceDepth::initSolver] Error, param vectors must be of same length.");

  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FitSurfaceDepth::initSolver] Error, surface not initialized (initSurface).");

  m_K = SparseMatrix( param0.rows(), m_nurbs.CVCount() );

  typedef Eigen::Triplet<double> Tri;
  std::vector<Tri> tripletList;
  tripletList.resize( param0.rows() * m_nurbs.Order(0) * m_nurbs.Order(1) );

  if(!m_quiet)
    printf("[FitSurfaceDepth::initSolver] entries: %lu  rows: %lu  cols: %d\n",
           param0.rows()* m_nurbs.Order(0) * m_nurbs.Order(1),
           param0.rows(),
           m_nurbs.CVCount());

  double *N0 = new double[m_nurbs.Order (0) * m_nurbs.Order (0)];
  double *N1 = new double[m_nurbs.Order (1) * m_nurbs.Order (1)];
  int E,F;
  int ti(0);

  for(Eigen::MatrixXd::Index row=0; row<param0.rows(); row++)
  {
    E = ON_NurbsSpanIndex (m_nurbs.m_order[0], m_nurbs.m_cv_count[0], m_nurbs.m_knot[0], param0(row), 0, 0);
    F = ON_NurbsSpanIndex (m_nurbs.m_order[1], m_nurbs.m_cv_count[1], m_nurbs.m_knot[1], param1(row), 0, 0);

    ON_EvaluateNurbsBasis (m_nurbs.Order (0), m_nurbs.m_knot[0] + E, param0(row), N0);
    ON_EvaluateNurbsBasis (m_nurbs.Order (1), m_nurbs.m_knot[1] + F, param1(row), N1);

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
    throw std::runtime_error("[FitSurfaceDepth::initSolver] decomposition failed.");

  if(!m_quiet)
    printf("[FitSurfaceDepth::initSolver] decomposition done\n");

  m_use_indices = false;
}

void FitSurfaceDepth::initSolver(const Eigen::VectorXd& param0,
                                 const Eigen::VectorXd& param1,
                                 const std::vector<int>& indices)
{
  if(param0.rows()!=param1.rows())
    throw std::runtime_error("[FitSurfaceDepth::initSolver] Error, param vectors must be of same length.");

  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FitSurfaceDepth::initSolver] Error, surface not initialized (initSurface).");

  m_K = SparseMatrix( indices.size(), m_nurbs.CVCount() );

  typedef Eigen::Triplet<double> Tri;
  std::vector<Tri> tripletList;
  tripletList.resize( indices.size() * m_nurbs.Order(0) * m_nurbs.Order(1) );

  if(!m_quiet)
    printf("[FitSurfaceDepth::initSolver] entries: %lu  rows: %lu  cols: %d\n",
           indices.size() * m_nurbs.Order(0) * m_nurbs.Order(1),
           indices.size(),
           m_nurbs.CVCount());

  double *N0 = new double[m_nurbs.Order (0) * m_nurbs.Order (0)];
  double *N1 = new double[m_nurbs.Order (1) * m_nurbs.Order (1)];
  int E,F;
  int ti(0);

  for(size_t row=0; row<indices.size(); row++)
  {
    const double& p0 = param0(indices[row]);
    const double& p1 = param1(indices[row]);

    E = ON_NurbsSpanIndex (m_nurbs.m_order[0], m_nurbs.m_cv_count[0], m_nurbs.m_knot[0], p0, 0, 0);
    F = ON_NurbsSpanIndex (m_nurbs.m_order[1], m_nurbs.m_cv_count[1], m_nurbs.m_knot[1], p1, 0, 0);

    ON_EvaluateNurbsBasis (m_nurbs.Order (0), m_nurbs.m_knot[0] + E, p0, N0);
    ON_EvaluateNurbsBasis (m_nurbs.Order (1), m_nurbs.m_knot[1] + F, p1, N1);

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
    throw std::runtime_error("[FitSurfaceDepth::initSolver] decomposition failed.");
  if(!m_quiet)
    printf("[FitSurfaceDepth::initSolver] decomposition done\n");
  m_use_indices = true;
}

void FitSurfaceDepth::solve(const Eigen::VectorXd &z)
{
  if(m_use_indices)
    throw std::runtime_error("[FitSurfaceDepth::solve] Error, solver initialized with indices (use solve(Eigen::VectorXd&, const std::vector<int>&) instead.\n");

  m_b = m_solver->solve(z);

  updateSurf();
}

void FitSurfaceDepth::solve(const Eigen::VectorXd &z, const std::vector<int>& indices)
{
  if(!m_use_indices)
    throw std::runtime_error("[FitSurfaceDepth::solve] Error, solver initialized without indices (use solve(Eigen::VectorXd&) instead.\n");

  Eigen::VectorXd zn(m_solver->rows(),1);

  for(size_t i=0; i<indices.size(); i++)
    zn(i) = z(indices[i]);

  m_b = m_solver->solve(zn);

  updateSurf();
}

void FitSurfaceDepth::updateSurf()
{
  int ncp = m_nurbs.CVCount ();

  if(m_b.rows()!=ncp)
    throw std::runtime_error("[FitSurfaceDepth::updateSurf] Error, number of control points does not match.");

  for (int i = 0; i < ncp; i++)
  {
    ON_3dPoint cp;

    cp.x = m_b(i, 0);
    cp.y = cp.x;
    cp.z = cp.x;

    m_nurbs.SetCV (gl2gr(i), gl2gc(i), cp);
  }
}

Eigen::VectorXd FitSurfaceDepth::GetError(const Eigen::VectorXd& z)
{
  if(m_use_indices)
    throw std::runtime_error("[FitSurfaceDepth::GetError] Error, solver initialized with indices (use solve(Eigen::VectorXd&, const std::vector<int>&) instead.\n");

  // compute K*b (i.e. points on surface)
  Eigen::VectorXd s(z.rows(),1);
  s.setZero();
  for (int k=0; k<m_K.outerSize(); ++k)
  for (SparseMatrix::InnerIterator it(m_K,k); it; ++it)
    s(it.row()) += it.value() * m_b(it.col());

  // return (K*b-z)
  return (s-z);
}

Eigen::VectorXd FitSurfaceDepth::GetError(const Eigen::VectorXd& z, const std::vector<int>& indices)
{
  if(!m_use_indices)
    throw std::runtime_error("[FitSurfaceDepth::GetError] Error, solver initialized without indices (use solve(Eigen::VectorXd&) instead.\n");

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
