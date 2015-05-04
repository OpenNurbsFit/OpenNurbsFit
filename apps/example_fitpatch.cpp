
#include "triangulation.h"

#include "TomGine/tgTomGineThread.h"
#include "NurbsFit/patch.h"

using namespace TomGine;
using namespace nurbsfit;

TomGine::tgTomGineThread viewer(800,600, "Example FitSurfaceDepth");

void CreateTent(Eigen::VectorXd& param0, Eigen::VectorXd& param1, Eigen::VectorXd& values,
                size_t num, FitSurface::Domain domain, int dim)
{
  param0 = Eigen::VectorXd(num);
  param1 = Eigen::VectorXd(num);
  values = Eigen::VectorXd(num*dim);

  for(size_t i=0; i<num; i++)
  {
    param0(i) = domain.width  * double(rand())/RAND_MAX;
    param1(i) = domain.height * double(rand())/RAND_MAX;

    if(dim==1)
    {
      if( param0(i) < domain.height*0.5 )
        values(i) = param0(i)/domain.height;
      else
        values(i) = 1.0 - param0(i)/domain.height;
    }

    if(dim==3)
    {
      values(dim*i+0) = param0(i) + domain.x;
      values(dim*i+1) = param1(i) + domain.y;

      if( param0(i) < domain.height*0.5 )
        values(dim*i+2) = param0(i)/domain.height;
      else
        values(dim*i+2) = 1.0 - param0(i)/domain.height;
    }

    param0(i) += domain.x;
    param1(i) += domain.y;
  }
}

void Nurbs2Mesh(ON_NurbsSurface& nurbs, TomGine::tgModel& mesh, int resX=32, int resY=32)
{
  double x0 = nurbs.Knot (0, 0);
  double x1 = nurbs.Knot (0, nurbs.m_knot_capacity[0] - 1);
  double y0 = nurbs.Knot (1, 0);
  double y1 = nurbs.Knot (1, nurbs.m_knot_capacity[1] - 1);

  mesh.Clear();
  tgShapeCreator::CreatePlaneXY(mesh, x0,y0,0, x1-x0,y1-y0, resX-1, resY-1);

  double z[nurbs.Dimension()];
  for(size_t i=0; i<mesh.m_vertices.size(); i++)
  {
    TomGine::tgVertex& v = mesh.m_vertices[i];
    nurbs.Evaluate (v.pos.x, v.pos.y, 0, nurbs.Dimension(), z);

    if(nurbs.Dimension()==3)
    {
      v.pos.x = z[0];
      v.pos.y = z[1];
      v.pos.z = z[2];
    }

    if(nurbs.Dimension()==1)
    {
      v.pos.z = z[0];
    }

  }

  mesh.ComputeNormals();
}

int main(int argc, char *argv[])
{
  viewer.SetClearColor(0.5f);

  // define domain
  FitSurface::Domain domain(0.0,0.0,1.0,1.0);

  // ###################### DIMENSION = 1 ######################
  // create data points
  int dim(3);
  Eigen::VectorXd param0, param1, values;
  srand(0);
  CreateTent(param0,param1,values, 1000, domain, dim);

  for(int i=0; i<param0.rows(); i++)
    viewer.AddPoint3D(values(3*i+0),values(3*i+1),values(3*i+2));

  // fit nurbs surface
  FitPatch fit;
  fit.initSurface(dim, 3,3, 10,10, domain);
  fit.initSolver(param0,param1);
  fit.solve(values);

  // visualize
  ON_NurbsSurface nurbs = fit.getSurface();
  TomGine::tgRenderModel mesh;
  mesh.SetColor(50,50,150);
  Triangulation::convertNurbs2tgModel(nurbs,mesh, 64,64);
  viewer.AddModel3D(mesh);

  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Space);

  for(size_t j=0; j<100; j++)
  {
    Eigen::Vector2d param;
    Eigen::VectorXd value(dim);
    int iter(0), k;
    double accuracy(0.0), a;
    for(unsigned i=0; i<param0.rows(); i++)
    {
      value(0) = values(dim*i + 0);
      value(1) = values(dim*i + 1);
      value(2) = values(dim*i + 2);
      param = fit.reparameterize(value, Eigen::Vector2d(param0(i),param1(i)), k, a, 100, 1e-6);
      iter += k;
      accuracy += a;
      param0(i) = param(0);
      param1(i) = param(1);
    }
    printf("steps: %d  accuracy: %f\n", iter, accuracy);
    fit.initSolver(param0,param1);
    fit.solve(values);

    if(accuracy<1e-4)
      break;

  }

  // visualize
  nurbs = fit.getSurface();
  mesh.SetColor(150,50,50);
  Triangulation::convertNurbs2tgModel(nurbs,mesh, 64,64);
  viewer.AddModel3D(mesh);


  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return 0;
}
