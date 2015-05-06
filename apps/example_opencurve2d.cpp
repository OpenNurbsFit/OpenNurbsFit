
#include "TomGine/tgTomGineThread.h"
#include "NurbsFit/open_curve.h"

using namespace TomGine;
using namespace nurbsfit;

TomGine::tgTomGineThread viewer(800,600, "Example: Fit Open Curve 2D");

// create data points for curve fitting
void CreateSine(Eigen::VectorXd& params, Eigen::VectorXd& points,
                int dims, int num_points, double k=1.0, double a=1.0, double b=1.0)
{
  params.resize(num_points);
  points.resize(num_points*dims);

  for(int i=0; i<num_points; i++)
  {
    double x = double(i) / (num_points - 1);
    double y = x * 2.0 * M_PI;
    params(i) = x;

    if(dims==1)
    {
      points(i) = a * sin(k*y);
    }

    if(dims==2)
    {
      points(i*dims) = x;
      points(i*dims+1) = a * sin(k*y);
    }

    if(dims==3)
    {
      points(i*dims) = x;
      points(i*dims+1) = a * sin(k*y);
      points(i*dims+2) = b * cos(k*y);
    }
  }
}

// Converts a B-spline curve to a mesh model
void Nurbs2Mesh(ON_NurbsCurve& nurbs, TomGine::tgModel& mesh, int res=256)
{
  double x0 = nurbs.Knot (0);
  double x1 = nurbs.Knot (nurbs.m_knot_capacity - 1);
  double w = x1-x0;

  mesh.Clear();

  mesh.m_lines.resize(res);

  double z0[nurbs.Dimension()];
  double z1[nurbs.Dimension()];
  for(int i=0; i<res-1; i++)
  {
    tgLine& line = mesh.m_lines[i];
    double param0 = x0 + w * double(i) / (res-1);
    double param1 = x0 + w * double(i+1) / (res-1);

    nurbs.Evaluate(param0, 0, nurbs.Dimension(), z0);
    nurbs.Evaluate(param1, 0, nurbs.Dimension(), z1);

    if(nurbs.Dimension()==3)
    {
      line.start.x = z0[0]; line.end.x = z1[0];
      line.start.y = z0[1]; line.end.y = z1[1];
      line.start.z = z0[2]; line.end.z = z1[2];
    }

    if(nurbs.Dimension()==2)
    {
      line.start.x = z0[0]; line.end.x = z1[0];
      line.start.y = z0[1]; line.end.y = z1[1];
      line.start.z = 0.0;   line.end.z = 0.0;
    }

    if(nurbs.Dimension()==1)
    {
      line.start.x = param0;  line.end.x = param1;
      line.start.y = z0[0];   line.end.y = z1[0];
      line.start.z = 0.0;     line.end.z = 0.0;
    }

  }

  if(nurbs.Dimension()==1)
  {
    double gx[nurbs.CVCount()];
    nurbs.GetGrevilleAbcissae(gx);

    for(int i=0; i<nurbs.CVCount(); i++)
    {
      ON_3dPoint cp;
      nurbs.GetCV(i,cp);
      mesh.m_points.push_back(vec3(gx[i],cp.x,0.0));
    }

  }else{

    for(int i=0; i<nurbs.CVCount(); i++)
    {
      ON_3dPoint cp;
      nurbs.GetCV(i,cp);
      mesh.m_points.push_back(vec3(cp.x,cp.y,cp.z));
    }
  }
  mesh.m_point_size = 10.0f;

}

int main(int argc, char *argv[])
{
  int num_points = 100; // number of data points
  int dims = 2;         // dimension of points & B-spline curve
  int cps = 10;         // number of control points
  int order = 4;        // polynomial order of B-spline curve
  bool clamp = false;   // clamp curve at ends, or leave them open

  Eigen::VectorXd params; // parameter values in B-spline domain, corresponding to points
  Eigen::VectorXd points; // data points B-spline is fitted to (params.size * dims == points.size)

  // create parameters and data points
  CreateSine(params, points, dims, num_points, 2.0, 0.25, 0.25);

  for(Eigen::VectorXd::Index i=0; i<num_points; i++)
  {
    if(dims==1)
      viewer.AddPoint3D(params(i), points(i), 0, 255,255,255, 5.0f);

    if(dims==2)
      viewer.AddPoint3D(points(i*dims), points(i*dims+1), 0, 255,255,255, 5.0f);

    if(dims==3)
      viewer.AddPoint3D(points(i*dims), points(i*dims+1), points(i*dims+2), 255,255,255, 5.0f);
  }

  // fit an open B-spline curve
  FitOpenCurve::Domain range(0,1);
  FitOpenCurve fit;
  fit.setQuiet(false);
  fit.initCurve(dims, order, cps, range, clamp);
  fit.initSolver(params);
  fit.solve(points);

  // visualize
  ON_NurbsCurve nurbs = fit.getCurve();

  ON_TextLog out;
  nurbs.Dump(out);

  tgModel mesh;
  Nurbs2Mesh(nurbs, mesh);
  mesh.m_line_width = 3.0f;
  mesh.m_point_color = vec3(0,1,0);
  viewer.AddModel3D(mesh);


  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return 0;
}
