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

#include "patch.h"
#include <stdexcept>

using namespace nurbsfit;

void FitPatch::initSurface(int dims, int order0, int order1, int cps0, int cps1, Domain roi)
{
  if( cps0 < order0 )
    cps0 = order0;

  if( cps1 < order1 )
    cps1 = order1;

  m_nurbs = ON_NurbsSurface (dims, false, order0, order1, cps0, cps1);

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

  m_x = Eigen::VectorXd(m_nurbs.CVCount()*dims,1);
  m_x.setZero();
  for (int i = 0; i < m_nurbs.CVCount(0); i++)
  {
    for (int j = 0; j < m_nurbs.CVCount(1); j++)
    {
      m_nurbs.SetCV (i, j, ON_3dPoint(0,0,0));
    }
  }
}

void FitPatch::initSolver(const Eigen::VectorXd& param0,
                                 const Eigen::VectorXd& param1)
{
  if(param0.rows()!=param1.rows())
    throw std::runtime_error("[FitPatch::initSolver] Error, param vectors must be of same length.");

  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FitPatch::initSolver] Error, surface not initialized (initSurface).");

  int dim = m_nurbs.Dimension();
  m_A = SparseMatrix( param0.rows()*dim, m_nurbs.CVCount()*dim );

  typedef Eigen::Triplet<double> Tri;
  std::vector<Tri> tripletList;
  tripletList.resize( 3*param0.rows() * m_nurbs.Order(0) * m_nurbs.Order(1) );

  if(!m_quiet)
    printf("[FitPatch::initSolver] entries: %lu  rows: %lu  cols: %d\n",
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
        tripletList[ti] = Tri( dim*row, dim*lrc2gl(E, F, i, j), N0[i] * N1[j] );
        if(dim>1)
          tripletList[ti+1] = Tri( dim*row+1, dim*lrc2gl(E, F, i, j)+1, N0[i] * N1[j] );
        if(dim>2)
          tripletList[ti+2] = Tri( dim*row+2, dim*lrc2gl(E, F, i, j)+2, N0[i] * N1[j] );
        ti+=dim;
      } // j
    } // i

  } // row

  delete [] N1;
  delete [] N0;

  m_A.setFromTriplets(tripletList.begin(), tripletList.end());

  if(m_solver==NULL)
    m_solver = new SPQR();
  m_solver->compute(m_A);
  if(m_solver->info()!=Eigen::Success)
    throw std::runtime_error("[FitPatch::initSolver] decomposition failed.");

  if(!m_quiet)
    printf("[FitPatch::initSolver] decomposition done\n");
}

void FitPatch::solve(const Eigen::VectorXd& values)
{
  m_x = m_solver->solve(values);

  updateSurf();
}

void FitPatch::updateSurf()
{
  int ncp = m_nurbs.CVCount ();

  if(m_x.rows()!=ncp*m_nurbs.Dimension())
    throw std::runtime_error("[FitPatch::updateSurf] Error, number of control points does not match.");

  for (int i = 0; i < ncp; i++)
  {
    ON_3dPoint cp;

    if(m_nurbs.Dimension()==1)
    {
      cp.x = m_x(i);
    }

    if(m_nurbs.Dimension()==2)
    {
      cp.x = m_x(2*i+0);
      cp.y = m_x(2*i+1);
    }

    if(m_nurbs.Dimension()==3)
    {
      cp.x = m_x(3*i+0);
      cp.y = m_x(3*i+1);
      cp.z = m_x(3*i+2);
    }

    m_nurbs.SetCV(gl2gr(i), gl2gc(i), cp);
  }
}

Eigen::VectorXd FitPatch::getError(const Eigen::VectorXd& values)
{
  // compute A*x (i.e. points on surface)
  Eigen::VectorXd Ax(values.rows(),1);
  Ax.setZero();
  for (int k=0; k<m_A.outerSize(); ++k)
  for (SparseMatrix::InnerIterator it(m_A,k); it; ++it)
    Ax(it.row()) += it.value() * m_x(it.col());

  // return (A*x-b)
  return (Ax-values);
}
