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
 */

#ifndef NURBS_FIT_CURVE_H
#define NURBS_FIT_CURVE_H

#include "OpenNurbs/opennurbs.h"

#undef Success
#include <Eigen/SparseQR>
#include <Eigen/SPQRSupport>

/** @brief nurbsfit namespace */
namespace nurbsfit{

/** @brief Fitting a NURBS curve to data. */
class FitCurve
{
public:
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::SPQR<SparseMatrix> SPQR;

  /** @brief Defines linear parametric domain of ON_NurbsCurve. */
  struct Domain
  {
    double start;    /// minimal point of domain
    double end;      /// maximal point of domain

    Domain() : start(0), end(1) {}

    /// Direct construction.
    Domain(double _start, double _end)
      : start(_start), end(_end) {}

    /// Construction of domain by bounding box of parametric positions of values (b).
    Domain(const Eigen::VectorXd &param);

    double width(){ return (end-start); }
  };

protected:
  bool m_quiet;             /// suppress console output
  ON_NurbsCurve m_nurbs;    /// curve fitted to data
  Eigen::VectorXd m_x;      /// control points
  SparseMatrix m_A;         /// matrix of linear system (containing basis functions of curve)
  SPQR* m_solver;           /// solving for x:  A * x = b

public:
  FitCurve() : m_quiet(true), m_solver(NULL) {}
  ~FitCurve(){ if(m_solver!=NULL) delete m_solver; }

  /// Accessor to NURBS curve.
  inline const ON_NurbsCurve& getCurve() const { return m_nurbs; }

  inline void setQuiet(bool b) { m_quiet = b; }

};

} // namespace nurbsfit

#endif // NURBS_FIT_CURVE_H
