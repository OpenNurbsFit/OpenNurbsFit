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

#ifndef NURBS_FIT_OPEN_CURVE_H
#define NURBS_FIT_OPEN_CURVE_H

#include "curve.h"

namespace nurbsfit{

class FitOpenCurve : public FitCurve
{
protected:
  void updateCurve();

public:

  /** Initialize the NURBS curve.
   *  @param[in] dims Dimension of NURBS surface (ie. control points).
   *  @param[in] order Polynomial order
   *  @param[in] cps Number of control points
   *  @param[in] domain The domain range of the NURBS curve.
   *  @see FitCurve::Domain
   */
  void initCurve(int dims, int order, int cps, Domain range);

  /** Initialize the solver. Assembly of matrix A and QR decomposition. Requires NURBS curve to be initialized.
   *  @param[in] param Vector of parametric positions, corresponding to values
   */
  void initSolver(const Eigen::VectorXd& param);

  /** Solve the linear system A * x = b with respect to x und updates control points of ON_NurbsCurve m_nurbs.
   *  Requires NURBS curve and solver to be initialized.
   *  @param[in] values The values the NURBS curve is fitted to (conforms b in linear system)
   */
  void solve(const Eigen::VectorXd& values);

  /** Compute fitting error.
   *  @param[in] values The values the NURBS curve has been fitted to (conforms b in linear system)
   *  @return Fitting error for each value (return A*x-b)
   */
  Eigen::VectorXd getError(const Eigen::VectorXd& values);
};

} // namespace nurbsfit

#endif // NURBS_FIT_OPEN_CURVE_H
