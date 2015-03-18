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

#include <Eigen/SparseQR>
#include <Eigen/SPQRSupport>

namespace nurbsfit
{

/** \brief Fitting a 1D NurbsSurface to 1D-data */
class FitSurfaceDepth : public FitSurface
{
public:

  struct ROI
  {
    double x, y, width, height;
    ROI() : x(0), y(0), width(0), height(0) {}
    ROI(double _x, double _y, double _w, double _h) : x(_x), y(_y), width(_w), height(_h) {}
    ROI(const Eigen::MatrixXd &points);
  };

protected:
  ROI m_roi;

  typedef Eigen::SparseMatrix<double> SparseMatrix;

  Eigen::VectorXd m_b; // control points
  SparseMatrix m_K;
  Eigen::VectorXd m_rs;

  bool m_quiet;

  typedef Eigen::SPQR<SparseMatrix> SPQR;
  SPQR* m_solver;

  bool m_use_indices;

  void updateSurf();

public:
  FitSurfaceDepth() : m_quiet(true), m_solver(NULL) {}
  FitSurfaceDepth(int order,
                      int cps_width, int cps_height,
                      ROI img_roi,
                      const Eigen::MatrixXd& points);
  FitSurfaceDepth(int order,
                      int cps_width, int cps_height,
                      ROI img_roi,
                      const Eigen::MatrixXd& points,
                      const std::vector<int>& indices);
  virtual ~FitSurfaceDepth();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  inline ON_NurbsSurface& getSurface ()
  {
    return m_nurbs;
  }

  void initSurface(int order, int cps_width, int cps_height, ROI img_roi);

  virtual void initSolver(const Eigen::MatrixXd& points);
  virtual void initSolver(const Eigen::MatrixXd& points, const std::vector<int>& indices);

  virtual void solve(const Eigen::VectorXd& z);
  virtual void solve(const Eigen::VectorXd& z, const std::vector<int>& indices);

  Eigen::VectorXd GetError(const Eigen::VectorXd& z);
  Eigen::VectorXd GetError(const Eigen::VectorXd& z, const std::vector<int>& indices);

};

}

#endif
