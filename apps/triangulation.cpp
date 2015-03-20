/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011, Thomas Mörwald, Jonathan Balzer, Inc.
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
 *   * Neither the name of Thomas Mörwald or Jonathan Balzer nor the names of its
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
 * @author thomas.moerwald
 *
 */

#include "triangulation.h"
#include "TomGine/tgShapeCreator.h"
#include <stdexcept>

using namespace nurbsfit;

void Triangulation::reverse (ON_NurbsCurve &curve, bool z_negative)
{
  ON_3dPoint p0;
  curve.GetCV (0, p0);

  double z (0.0);

  for (int i = 1; i < curve.CVCount () - 1; i++)
  {
    ON_3dPoint p1, p2;
    curve.GetCV (i, p1);
    curve.GetCV (i + 1, p2);

    z += p1.x * p2.y - p1.y * p2.x;
  }

  if ((z_negative && z <= 0.0) || (!z_negative && z >= 0))
  {
    ON_NurbsCurve curve2 = curve;
    for (int i = 0; i < curve.CVCount (); i++)
    {
      int j = curve.CVCount () - 1 - i;
      ON_3dPoint p;
      curve.GetCV (i, p);
      curve2.SetCV (j, p);
    }
    curve = curve2;
  }
}

void Triangulation::flip (int dir, ON_NurbsSurface &nurbs)
{
  ON_NurbsSurface result (nurbs);

  for (int i = 0; i < nurbs.CVCount (0); i++)
  {
    for (int j = 0; j < nurbs.CVCount (1); j++)
    {
      ON_3dPoint cp;
      nurbs.GetCV (i, j, cp);
      if (dir == 0)
        result.SetCV (nurbs.CVCount (0) - i - 1, j, cp);
      else if (dir == 1)
        result.SetCV (i, nurbs.CVCount (1) - j - 1, cp);
    }
  }

  nurbs = result;
}

void Triangulation::convertNurbs2tgModel (const ON_NurbsSurface &nurbs, TomGine::tgModel &model,
                                          unsigned resU, unsigned resV, bool cps)
{
  if (nurbs.m_knot_capacity[0] <= 1 || nurbs.m_knot_capacity[1] <= 1)
    throw std::runtime_error("[Triangulation::convertNurbs2tgModel] Warning: ON knot vector empty.\n");

  model.Clear ();

  double x0 = nurbs.Knot (0, 0);
  double x1 = nurbs.Knot (0, nurbs.m_knot_capacity[0] - 1);
  double w = x1 - x0;
  double y0 = nurbs.Knot (1, 0);
  double y1 = nurbs.Knot (1, nurbs.m_knot_capacity[1] - 1);
  double h = y1 - y0;

  TomGine::tgShapeCreator::CreatePlaneXY (model, x0, y0, 0.0, w, h, resU, resV);
  TomGine::vec3 tu, tv;
  double pointAndTangents[9];

  for (unsigned i = 0; i < model.m_vertices.size (); i++)
  {

    TomGine::tgVertex& v = model.m_vertices[i];

    if(nurbs.Dimension()==3)
    {
      nurbs.Evaluate (v.pos.x, v.pos.y, 1, nurbs.Dimension(), pointAndTangents);

      v.pos.x = pointAndTangents[0];
      v.pos.y = pointAndTangents[1];
      v.pos.z = pointAndTangents[2];

      tu.x = pointAndTangents[3]; // use tu
      tu.y = pointAndTangents[4];
      tu.z = pointAndTangents[5];
      tv.x = pointAndTangents[6]; // use tv
      tv.y = pointAndTangents[7];
      tv.z = pointAndTangents[8];

      v.normal = TomGine::normalize(TomGine::cross (tu, tv));
    }

    if(nurbs.Dimension()==1)
    {
      nurbs.Evaluate (v.pos.x, v.pos.y, 0, nurbs.Dimension(), &pointAndTangents[0]);
      v.pos.z = pointAndTangents[0];
    }
  }

  if(nurbs.Dimension()==1)
    model.ComputeNormals();


  if (cps)
  {
    model.m_points.clear ();
    model.m_point_size = 2.0f;
    for (int j = 0; j < nurbs.CVCount (1); j++)
    {
      for (int i = 0; i < nurbs.CVCount (0); i++)
      {
        ON_3dPoint p;
        nurbs.GetCV (i, j, p);
        model.m_points.push_back (TomGine::vec3 (p.x, p.y, p.z));
      }
    }
  }
}

