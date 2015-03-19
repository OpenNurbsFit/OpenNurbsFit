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

#ifndef NURBS_FIT_SURFACE_DEPTH_H
#define NURBS_FIT_SURFACE_DEPTH_H

#include "surface.h"

namespace nurbsfit
{

/** \brief Fitting a 1D NurbsSurface to 1D-data */
class FitPatch : public FitSurface
{
protected:
  void updateSurf();

  // global row/col index to global lexicographic index
  int grc2gl (int I, int J) const { return m_nurbs.CVCount (1) * I + J; }
  // local row/col index to global lexicographic index
  int lrc2gl (int E, int F, int i, int j) const { return grc2gl (E + i, F + j); }
  // global lexicographic in global row index
  int gl2gr (int A) const { return (static_cast<int> (A / m_nurbs.CVCount (1))); }
  // global lexicographic in global col index
  int gl2gc (int A) const { return (static_cast<int> (A % m_nurbs.CVCount (1))); }

public:

  void initSurface(int dims, int order0, int order1, int cps0, int cps1, Domain roi);

  void initSolver(const Eigen::VectorXd& param0, const Eigen::VectorXd& param1);

  void solve(const Eigen::VectorXd& values);

//  Eigen::VectorXd GetError(const Eigen::VectorXd& z);
//  Eigen::VectorXd GetError(const Eigen::VectorXd& z, const std::vector<int>& indices);

};

}

#endif
