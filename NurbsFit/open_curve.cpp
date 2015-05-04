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

#include "open_curve.h"
#include <stdexcept>

using namespace nurbsfit;

void FitOpenCurve::initCurve(int dims, int order, int cps, Domain range)
{
  if(cps<order)
    cps = order;

  m_nurbs = ON_NurbsCurve(dims, false, order, cps);

  double dd = range.width() / (m_nurbs.KnotCount() - 2*(order-2) - 1);

  m_nurbs.MakeClampedUniformKnotVector(dd);

  for (int i = 0; i < m_nurbs.KnotCount(); i++)
    m_nurbs.SetKnot (i, m_nurbs.Knot(i) + range.start);

  m_x = Eigen::VectorXd(m_nurbs.CVCount()*dims,1);
  m_x.setZero();
  for (int i = 0; i < m_nurbs.CVCount(); i++)
    m_nurbs.SetCV(i, ON_3dPoint(m_nurbs.Knot(i),0,0));
}

void FitOpenCurve::initSolver(const Eigen::VectorXd& params)
{
  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FitOpenCurve::initSolver] Error, curve not initialized (initCurve).");

  if(params.rows()<2)
    throw std::runtime_error("[FitOpenCurve::initSolver] Error, insufficient parameter points (<2).");

  int dim = m_nurbs.Dimension();
  m_A = SparseMatrix( params.rows()*dim, m_nurbs.CVCount()*dim );

  typedef Eigen::Triplet<double> Tri;
  std::vector<Tri> tripletList;
  tripletList.resize( dim * params.rows() * m_nurbs.Order() );

  if(!m_quiet)
    printf("[FitOpenCurve::initSolver] entries: %lu  rows: %lu  cols: %d\n",
           dim * params.rows()* m_nurbs.Order(),
           dim * params.rows(),
           m_nurbs.CVCount());

  double *N = new double[m_nurbs.m_order * m_nurbs.m_order];
  int E;
  int ti(0);

  for(Eigen::MatrixXd::Index row=0; row<params.rows(); row++)
  {
    E = ON_NurbsSpanIndex (m_nurbs.m_order, m_nurbs.m_cv_count, m_nurbs.m_knot, params(row), 0, 0);

    ON_EvaluateNurbsBasis (m_nurbs.m_order, m_nurbs.m_knot + E, params(row), N);

    for (int i = 0; i < m_nurbs.Order(); i++)
    {
        tripletList[ti] = Tri( dim*row, dim*(E+i), N[i] );
        if(dim>1)
          tripletList[ti+1] = Tri( dim*row+1, dim*(E+i)+1, N[i] );
        if(dim>2)
          tripletList[ti+2] = Tri( dim*row+2, dim*(E+i)+2, N[i] );

        ti+=dim;
    } // i

  } // row

  delete [] N;

  m_A.setFromTriplets(tripletList.begin(), tripletList.end());

  if(m_solver==NULL)
    m_solver = new SPQR();
  m_solver->compute(m_A);
  if(m_solver->info()!=Eigen::Success)
    throw std::runtime_error("[FitOpenCurve::initSolver] decomposition failed.");

  if(!m_quiet)
    printf("[FitOpenCurve::initSolver] decomposition done\n");
}

void FitOpenCurve::solve(const Eigen::VectorXd& values)
{
  m_x = m_solver->solve(values);

  updateCurve();
}

void FitOpenCurve::updateCurve()
{
  int ncp = m_nurbs.CVCount();

  if(m_x.rows()!=ncp*m_nurbs.Dimension())
    throw std::runtime_error("[FitOpenCurve::updateSurf] Error, number of control points does not match.");

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

    m_nurbs.SetCV(i, cp);
  }
}
