
#include "TomGine/tgTomGineThread.h"
#include "NurbsFit/open_curve.h"

using namespace TomGine;
using namespace nurbsfit;

TomGine::tgTomGineThread viewer(800,600, "Example: Fit Open Curve 2D");

// create data points for curve fitting
void CreateSine(Eigen::VectorXd& params, Eigen::VectorXd& points,
                FitOpenCurve::Domain& range,
                int dims, int num_points, double k=1.0, double a=1.0, double b=1.0)
{
  params.resize(num_points);
  points.resize(num_points*dims);

  for(int i=0; i<num_points; i++)
  {
    double x = double(i) / (num_points - 1);
    double y = x * 2.0 * M_PI;
    params(i) = range.width() * x;

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
void Nurbs2MeshManually(ON_NurbsCurve& nurbs, TomGine::tgModel& mesh, int res=256,
                        FitOpenCurve::Domain range = FitOpenCurve::Domain(0,0), double xoff=0.001)
{
  if(range.start==range.end)
  {
    range.start = nurbs.Knot(0);
    range.end = nurbs.Knot(nurbs.KnotCount()-1);
  }

  double w = range.width();

  mesh.Clear();

  mesh.m_lines.resize(res);

  double tangent_scale = 1.0;
  int nder = 1;
  int nvals = nurbs.m_dim*(nder+1);
  double z0[nvals];
  double z1[nvals];

  int E0, E1;
  double *N0 = new double[nurbs.m_order * nurbs.m_order];
  double *N1 = new double[nurbs.m_order * nurbs.m_order];

  std::vector<ON_3dPoint> cps(nurbs.CVCount());
  for(int i=0; i<nurbs.CVCount(); i++)
    nurbs.GetCV(i, cps[i]);

  for(int i=0; i<res-1; i++)
  {
    tgLine& line = mesh.m_lines[i];
    double param0 = range.start + w * double(i) / (res-1);
    double param1 = range.start + w * double(i+1) / (res-1);

    if(nurbs.Dimension()==3)
    {
      // test: try to manually evaluate curve + derivatives
      E0 = ON_NurbsSpanIndex (nurbs.m_order, nurbs.m_cv_count, nurbs.m_knot, param0, 0, 0);
      ON_EvaluateNurbsBasis (nurbs.m_order, nurbs.m_knot + E0, param0, N0);
      ON_EvaluateNurbsBasisDerivatives(nurbs.m_order, nurbs.m_knot + E0, nder, N0);

      E1 = ON_NurbsSpanIndex (nurbs.m_order, nurbs.m_cv_count, nurbs.m_knot, param1, 0, 0);
      ON_EvaluateNurbsBasis (nurbs.m_order, nurbs.m_knot + E1, param1, N1);
      ON_EvaluateNurbsBasisDerivatives(nurbs.m_order, nurbs.m_knot + E1, nder, N1);

      for(int j=0; j<nvals; j++)
      {
        z0[j] = 0.0;
        z1[j] = 0.0;
      }

      for (int j = 0; j < nurbs.m_order; j++)
      {
        z0[0] += N0[j] * cps[E0+j].x;
        z0[1] += N0[j] * cps[E0+j].y;
        z0[2] += N0[j] * cps[E0+j].z;
        z1[0] += N1[j] * cps[E1+j].x;
        z1[1] += N1[j] * cps[E1+j].y;
        z1[2] += N1[j] * cps[E1+j].z;

        z0[3] += N0[nurbs.m_order+j] * cps[E0+j].x;
        z0[4] += N0[nurbs.m_order+j] * cps[E0+j].y;
        z0[5] += N0[nurbs.m_order+j] * cps[E0+j].z;
        z1[3] += N1[nurbs.m_order+j] * cps[E1+j].x;
        z1[4] += N1[nurbs.m_order+j] * cps[E1+j].y;
        z1[5] += N1[nurbs.m_order+j] * cps[E1+j].z;
      }

      line.start.x = xoff+z0[0]; line.end.x = xoff+z1[0];
      line.start.y = z0[1]; line.end.y = z1[1];
      line.start.z = z0[2]; line.end.z = z1[2];

      tgLine tangent;
      tangent.start.x = xoff+z0[0]; tangent.end.x = xoff+z0[0] + z0[3] * tangent_scale;
      tangent.start.y = z0[1]; tangent.end.y = z0[1] + z0[4] * tangent_scale;
      tangent.start.z = z0[2]; tangent.end.z = z0[2] + z0[5] * tangent_scale;
      mesh.m_lines.push_back(tangent);

    }
  }

  delete N0;
  delete N1;
}

// Converts a B-spline curve to a mesh model
void Nurbs2Mesh(ON_NurbsCurve& nurbs, TomGine::tgModel& mesh, int res=256,
                FitOpenCurve::Domain range = FitOpenCurve::Domain(0,0))
{
  if(range.start==range.end)
  {
    range.start = nurbs.Knot(0);
    range.end = nurbs.Knot(nurbs.KnotCount()-1);
  }

  double w = range.width();

  mesh.Clear();

  mesh.m_lines.resize(res);

  double tangent_scale = 1.0;
  int nder = 1;
  int nvals = nurbs.m_dim*(nder+1);
  double z0[nvals];
  double z1[nvals];

  for(int i=0; i<res-1; i++)
  {
    tgLine& line = mesh.m_lines[i];
    double param0 = range.start + w * double(i) / (res-1);
    double param1 = range.start + w * double(i+1) / (res-1);

    nurbs.Evaluate(param0, nder, nurbs.Dimension(), z0);
    nurbs.Evaluate(param1, nder, nurbs.Dimension(), z1);

    if(nurbs.Dimension()==3)
    {
      line.start.x = z0[0]; line.end.x = z1[0];
      line.start.y = z0[1]; line.end.y = z1[1];
      line.start.z = z0[2]; line.end.z = z1[2];
      {
        tgLine tangent;
        tangent.start.x = z0[0]; tangent.end.x = z0[0] + z0[3] * tangent_scale;
        tangent.start.y = z0[1]; tangent.end.y = z0[1] + z0[4] * tangent_scale;
        tangent.start.z = z0[2]; tangent.end.z = z0[2] + z0[5] * tangent_scale;
        mesh.m_lines.push_back(tangent);
      }

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
  int dims = 3;         // dimension of points & B-spline curve
  int cps = 15;         // number of control points
  int order = 4;        // polynomial order of B-spline curve
  bool clamp = false;   // clamp curve at ends, or leave them open

  Eigen::VectorXd params; // parameter values in B-spline domain, corresponding to points
  Eigen::VectorXd points; // data points B-spline is fitted to (params.size * dims == points.size)

  // create parameters and data points
  FitOpenCurve::Domain range(0,10);
  CreateSine(params, points, range, dims, num_points, 3.0, 0.25, 0.25);

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
  FitOpenCurve fit;
  fit.initCurve(dims, order, cps, range, clamp);
  fit.initSolver(params);
  fit.solve(points);

  // visualize
  ON_NurbsCurve nurbs = fit.getCurve();
  tgModel mesh;
  Nurbs2Mesh(nurbs, mesh,params.rows(),range);
  mesh.m_line_width = 3.0f;
  mesh.m_point_color = vec3(0,1,0);
  viewer.AddModel3D(mesh);

  Nurbs2MeshManually(nurbs, mesh,params.rows(),range);
  mesh.m_line_width = 3.0f;
  mesh.m_line_color = vec3(0.5,0.5,1.0);
  viewer.AddModel3D(mesh);


  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return 0;
}
