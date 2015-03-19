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
class FitSurfaceDepth : public FitSurface
{
protected:
  bool m_use_indices;

  void updateSurf();

public:
  FitSurfaceDepth(int order0, int order1,
                  int cps0, int cps1,
                  Domain roi,
                  const Eigen::MatrixXd& points);
  FitSurfaceDepth(int order0, int order1,
                  int cps0, int cps1,
                  Domain roi,
                  const Eigen::MatrixXd& points,
                  const std::vector<int>& indices);
  virtual ~FitSurfaceDepth();

  void initSurface(int order0, int order1, int cps0, int cps1, Domain roi);

  virtual void initSolver(const Eigen::VectorXd& param0,
                          const Eigen::VectorXd& param1);
  virtual void initSolver(const Eigen::VectorXd& param0,
                          const Eigen::VectorXd& param1,
                          const std::vector<int>& indices);

  virtual void solve(const Eigen::VectorXd& z);
  virtual void solve(const Eigen::VectorXd& z, const std::vector<int>& indices);

  Eigen::VectorXd GetError(const Eigen::VectorXd& z);
  Eigen::VectorXd GetError(const Eigen::VectorXd& z, const std::vector<int>& indices);

};

}

#endif
