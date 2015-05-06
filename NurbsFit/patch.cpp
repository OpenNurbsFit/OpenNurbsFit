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

#undef Success
#include <Eigen/Dense>

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
      m_nurbs.SetCV (i, j, ON_3dPoint(m_nurbs.Knot(0, i),m_nurbs.Knot(1, j),0));
    }
  }
}

void FitPatch::initSolver(const Eigen::VectorXd& param0, const Eigen::VectorXd& param1)
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

Eigen::Vector2d FitPatch::reparameterize(const Eigen::VectorXd &value, const Eigen::Vector2d &hint,
                                         int& steps, double &accuracy, int maxSteps, double minAccuracy)
{
  if(m_nurbs.CVCount() <= 0)
    throw std::runtime_error("[FitPatch::reparameterize] Error, surface not initialized (initSurface).");

  int nder = 1; // number of derivatives
  int dims = m_nurbs.Dimension();
  int nvals = dims*(nder+1)*(nder+2)/2;
  double pointAndTangents[nvals];

  Eigen::Vector2d current, delta, b, c;
  Eigen::Matrix2d A, I;
  Eigen::VectorXd p(dims), tu(dims), tv(dims), r(dims);
  I = Eigen::Matrix2d::Identity();

  double minU = m_nurbs.Knot(0,0);
  double minV = m_nurbs.Knot(1,0);
  double maxU = m_nurbs.Knot(0,m_nurbs.KnotCount(0)-1);
  double maxV = m_nurbs.Knot(1,m_nurbs.KnotCount(1)-1);
  double nu = 0.5;
  double lambda = 1.0;
  double accuracy_old(DBL_MAX);

  current = hint;

  for (steps = 0; steps < maxSteps; steps++)
  {

    m_nurbs.Evaluate (current(0), current(1), nder, dims, pointAndTangents);

    if(dims==1)
    {
      p(0) = pointAndTangents[0];
      tu(0) = pointAndTangents[1];
      tv(0) = pointAndTangents[2];
    }

    if(dims==2)
    {
      p(0) = pointAndTangents[0];
      p(1) = pointAndTangents[1];
      tu(0) = pointAndTangents[2];
      tu(1) = pointAndTangents[3];
      tv(0) = pointAndTangents[4];
      tv(1) = pointAndTangents[5];
    }

    if(dims==3)
    {
      p(0) = pointAndTangents[0];
      p(1) = pointAndTangents[1];
      p(2) = pointAndTangents[2];
      tu(0) = pointAndTangents[3];
      tu(1) = pointAndTangents[4];
      tu(2) = pointAndTangents[5];
      tv(0) = pointAndTangents[6];
      tv(1) = pointAndTangents[7];
      tv(2) = pointAndTangents[8];
    }

    r = p - value;

    b(0) = -r.dot (tu);
    b(1) = -r.dot (tv);

    A(0, 0) = tu.dot (tu);
    A(0, 1) = tu.dot (tv);
    A(1, 0) = A (0, 1);
    A(1, 1) = tv.dot (tv);

    delta = A.ldlt().solve(b);

    accuracy = delta.norm();
    if (accuracy < minAccuracy)
      return current;

    // step width control (quite heuristic)
    if(accuracy>accuracy_old*nu)
      lambda *= nu;
    accuracy_old = accuracy;

    // make step
    c = current + lambda * delta;

    // clamp to domain borders
    if (c(0) < minU)
      c(0) = minU;
    else if (c (0) > maxU)
      c(0) = maxU;
    if (c(1) < minV)
      c(1) = minV;
    else if (c(1) > maxV)
      c(1) = maxV;

    // compute real step
    delta = c - current;
    current = c;

    accuracy = delta.norm();
    if (accuracy < minAccuracy)
      return current;
  }

  printf ("[FitPatch::reparameterize] Warning: Method did not converge (%e %e %d)\n",
          accuracy, minAccuracy, maxSteps);
  printf ("[FitPatch::reparameterize]   (0: %f %f) (1: %f %f) %f %f ... %f %f\n",
          minU, maxU, minV, maxV,
          hint (0), hint (1), current (0), current (1));

  return current;
}

//int nvals = nurbs.m_dim*(nder+1)*(nder+2)/2;
//double P[nvals];
//nurbs.Evaluate (u, v, nder, nurbs.m_dim, P);

//// positions
//xx = P[0];    xy = P[1];    xz = P[2];

//// 1st derivatives (for normals)
//xu(0) = P[3];    xu(1) = P[4];    xu(2) = P[5];
//xv(0) = P[6];    xv(1) = P[7];    xv(2) = P[8];

//n = xu.cross(xv);
//n.normalize();
//v.normal = TomGine::vec3(n(0),n(1),n(2));

//// 2nd derivatives (for curvature)
//xuu(0) = P[9];     xuu(1) = P[10];    xuu(2) = P[11];
//xuv(0) = P[12];    xuv(1) = P[13];    xuv(2) = P[14];
//xvv(0) = P[15];    xvv(1) = P[16];    xvv(2) = P[17];

//#include <unsupported/Eigen/NumericalDiff>
//#include <unsupported/Eigen/NonLinearOptimization>

//struct SurfaceFunctor
//{
//  const Eigen::VectorXd& m_residuals;
//  const ON_NurbsSurface& m_nurbs;

//  SurfaceFunctor(const Eigen::VectorXd& residuals, const ON_NurbsSurface& nurbs)
//    : m_residuals(residuals), m_nurbs(nurbs){}

//  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& residuals) const
//  {


//    return 0;
//  }

//  int dx(const Eigen::VectorXd& x, Eigen::MatrixXd& jacobian)
//  {

//    return 0;
//  }

//  int inputs() const { return m_nurbs.CVCount(); }
//  int values() const { return m_residuals.rows(); }
//};

//Eigen::Vector2d FitPatch::reparameterizeLM(const Eigen::VectorXd &value, const Eigen::Vector2d &hint,
//                                         int& steps, double &accuracy, int maxSteps, double minAccuracy)
//{
//  if(m_nurbs.CVCount() <= 0)
//    throw std::runtime_error("[FitPatch::reparameterize] Error, surface not initialized (initSurface).");

//  int nder = 1; // number of derivatives
//  int dims = m_nurbs.Dimension();
//  int nvals = dims*(nder+1)*(nder+2)/2;
//  double pointAndTangents[nvals];

//  Eigen::Vector2d current, delta, b, c;
//  Eigen::Matrix2d A, I;
//  Eigen::VectorXd p(dims), tu(dims), tv(dims), r(dims);
//  I = Eigen::Matrix2d::Identity();

//  double minU = m_nurbs.Knot(0,0);
//  double minV = m_nurbs.Knot(1,0);
//  double maxU = m_nurbs.Knot(0,m_nurbs.KnotCount(0)-1);
//  double maxV = m_nurbs.Knot(1,m_nurbs.KnotCount(1)-1);

//  current = hint;

//  SurfaceFunctor func;
//  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor>,double> lm(numDiff);

//  for (steps = 0; steps < maxSteps; steps++)
//  {

//    m_nurbs.Evaluate (current(0), current(1), nder, dims, pointAndTangents);

//    if(dims==1)
//    {
//      p(0) = pointAndTangents[0];
//      tu(0) = pointAndTangents[1];
//      tv(0) = pointAndTangents[2];
//    }

//    if(dims==2)
//    {
//      p(0) = pointAndTangents[0];
//      p(1) = pointAndTangents[1];
//      tu(0) = pointAndTangents[2];
//      tu(1) = pointAndTangents[3];
//      tv(0) = pointAndTangents[4];
//      tv(1) = pointAndTangents[5];
//    }

//    if(dims==3)
//    {
//      p(0) = pointAndTangents[0];
//      p(1) = pointAndTangents[1];
//      p(2) = pointAndTangents[2];
//      tu(0) = pointAndTangents[3];
//      tu(1) = pointAndTangents[4];
//      tu(2) = pointAndTangents[5];
//      tv(0) = pointAndTangents[6];
//      tv(1) = pointAndTangents[7];
//      tv(2) = pointAndTangents[8];
//    }

//    r = p - value;

//    b(0) = -r.dot (tu);
//    b(1) = -r.dot (tv);

//    A(0, 0) = tu.dot (tu);
//    A(0, 1) = tu.dot (tv);
//    A(1, 0) = A (0, 1);
//    A(1, 1) = tv.dot (tv);

//    delta = A.ldlt().solve(b);

//    accuracy = delta.norm();
//    if (accuracy < minAccuracy)
//      return current;

//    // make step
//    c = current + delta;

//    // clamp to domain borders
//    if (c(0) < minU)
//      c(0) = minU;
//    else if (c (0) > maxU)
//      c(0) = maxU;
//    if (c(1) < minV)
//      c(1) = minV;
//    else if (c(1) > maxV)
//      c(1) = maxV;

//    // compute real step
//    delta = c - current;
//    current = c;

//    accuracy = delta.norm();
//    if (accuracy < minAccuracy)
//      return current;
//  }

//  printf ("[FitPatch::reparameterize] Warning: Method did not converge (%e %e %d)\n",
//          accuracy, minAccuracy, maxSteps);
//  printf ("[FitPatch::reparameterize]   (0: %f %f) (1: %f %f) %f %f ... %f %f\n",
//          minU, maxU, minV, maxV,
//          hint (0), hint (1), current (0), current (1));

//  return current;
//}

