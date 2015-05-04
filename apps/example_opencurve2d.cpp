
#include "TomGine/tgTomGineThread.h"
#include "NurbsFit/open_curve.h"

using namespace TomGine;
using namespace nurbsfit;

TomGine::tgTomGineThread viewer(800,600, "Example FitSurfaceDepth");

void CreateSine(Eigen::VectorXd& params, Eigen::VectorXd& values,
                int dims, int num_samples, double k=1.0, double a=1.0, double b=1.0)
{
  params.resize(num_samples);
  values.resize(num_samples*dims);

  for(int i=0; i<num_samples; i++)
  {
    double x = double(i) / (num_samples - 1);
    double y = x * 2.0 * M_PI;
    params(i) = x;

    if(dims==1)
    {
      values(i) = sin(y);
    }
    if(dims==2)
    {
      values(i*dims) = x;
      values(i*dims+1) = a * sin(k*y);
    }
    if(dims==3)
    {
      values(i*dims) = x;
      values(i*dims+1) = a * sin(k*y);
      values(i*dims+2) = b * cos(k*y);
    }
  }
}

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
      line.start.y = z0[1];   line.end.y = z1[1];
      line.start.z = 0.0;     line.end.z = 0.0;
    }

  }

  for(int i=0; i<nurbs.CVCount(); i++)
  {
    ON_3dPoint cp;
    nurbs.GetCV(i,cp);
    mesh.m_points.push_back(vec3(cp.x,cp.y,cp.z));
  }
  mesh.m_point_size = 10.0f;

}

int main(int argc, char *argv[])
{
  Eigen::VectorXd params;
  Eigen::VectorXd values;
  int num_samples = 1000;
  int dims = 2;
  int cps = 20;
  int order = 3;
  CreateSine(params, values, dims, num_samples, 2.0, 0.25, 0.25);

  for(Eigen::VectorXd::Index i=0; i<num_samples; i++)
  {
    if(dims==2)
      viewer.AddPoint3D(values(i*dims), values(i*dims+1), 0, 255,255,255, 5.0f);

    if(dims==3)
      viewer.AddPoint3D(values(i*dims), values(i*dims+1), values(i*dims+2), 255,255,255, 5.0f);
  }

  FitOpenCurve::Domain range(0,1);
  FitOpenCurve fit;
  fit.setQuiet(false);
  fit.initCurve(dims, order, cps, range);
  fit.initSolver(params);
  fit.solve(values);

  ON_NurbsCurve nurbs = fit.getCurve();
  tgModel mesh;
  Nurbs2Mesh(nurbs, mesh);
  mesh.m_line_width = 3.0f;

  viewer.AddModel3D(mesh);


  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return 0;
}
