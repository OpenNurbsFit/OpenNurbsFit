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

#ifndef NURBS_FIT_SURFACE_PATCH_H
#define NURBS_FIT_SURFACE_PATCH_H

#include "surface.h"

namespace nurbsfit
{

/** \brief Fitting a NURBS surface patch (clamped at borders) to data */
class FitPatch : public FitSurface
{
protected:
  void updateSurf();

  /// global row/col index to global lexicographic index
  int grc2gl (int I, int J) const { return m_nurbs.CVCount (1) * I + J; }
  /// local row/col index to global lexicographic index
  int lrc2gl (int E, int F, int i, int j) const { return grc2gl (E + i, F + j); }
  /// global lexicographic in global row index
  int gl2gr (int A) const { return (static_cast<int> (A / m_nurbs.CVCount (1))); }
  /// global lexicographic in global col index
  int gl2gc (int A) const { return (static_cast<int> (A % m_nurbs.CVCount (1))); }

public:

  /** Initialize the NURBS surface.
   *  @param[in] dims Dimension of NURBS surface (ie. control points).
   *  @param[in] order0 Polynomial order in x-direction.
   *  @param[in] order1 Polynomial order in y-direction.
   *  @param[in] cps0 Number of control points in x-direction.
   *  @param[in] cps1 Number of control poitns in y-direction.
   *  @param[in] domain The domain of the NURBS surface.
   *  @see FitSurface::Domain
   */
  void initSurface(int dims, int order0, int order1, int cps0, int cps1, Domain roi);

  /** Initialize the solver. Assembly of matrix K and QR decomposition. Requires NURBS surface to be initialized.
   *  @param[in] param0 Vector of parametric positions of values in x-dimension
   *  @param[in] param1 Vector of parametric positions of values in y-dimension
   */
  void initSolver(const Eigen::VectorXd& param0, const Eigen::VectorXd& param1);

  /** Solve the linear system A * x = b with respect to x und updates control points of ON_NurbsSurface m_nurbs.
   *  Requires NURBS surface and solver to be initialized.
   *  @param[in] values The values the NURBS surface is fitted to (conforms b in linear system)
   */
  void solve(const Eigen::VectorXd& values);

  /** Compute fitting error.
   *  @param[in] values The values the NURBS surface has been fitted to (conforms b in linear system)
   *  @return Fitting error for each value (return A*x-b)
   */
  Eigen::VectorXd getError(const Eigen::VectorXd& values);


  Eigen::Vector2d reparameterize(const Eigen::VectorXd &value, const Eigen::Vector2d &hint,
                                 int& steps, double& accuracy, int maxSteps, double minAccuracy);

  Eigen::Vector2d reparameterizeLM(const Eigen::VectorXd &value, const Eigen::Vector2d &hint,
                                   int& steps, double& accuracy, int maxSteps, double minAccuracy);

};

}

#endif
