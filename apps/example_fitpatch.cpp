
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

void Nurbs2Mesh(ON_NurbsSurface& nurbs, TomGine::tgModel& mesh,
                FitSurface::Domain domain, int resX=32, int resY=32)
{
  mesh.Clear();
  tgShapeCreator::CreatePlaneXY(mesh, domain.x,domain.y,0, domain.width,domain.height, resX-1, resY-1);

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
  {
    int dim(1);
    Eigen::VectorXd param0, param1, values;
    srand(0);
    CreateTent(param0,param1,values, 1000, domain, dim);

//    for(int i=0; i<param0.rows(); i++)
//      viewer.AddPoint3D(param0(i),param1(i),values(i,0));

    // fit nurbs surface
    FitPatch fit;
    fit.initSurface(dim, 3,3, 10,10, domain);
    fit.initSolver(param0,param1);
    fit.solve(values);

    // visualize
    ON_NurbsSurface nurbs = fit.getSurface();
    TomGine::tgRenderModel mesh;
    mesh.SetColor(200,50,50);
    Nurbs2Mesh(nurbs,mesh,domain);
    viewer.AddModel3D(mesh);
  }

  // ###################### DIMENSION = 3 ######################
  // create data points
  {
    int dim(3);
    Eigen::VectorXd param0, param1, values;
    srand(0);
    CreateTent(param0,param1,values, 1000, domain, dim);

//    for(int i=0; i<param0.rows(); i++)
//      viewer.AddPoint3D(param0(i),param1(i),values(i,0));

    // fit nurbs surface
    FitPatch fit;
    fit.initSurface(dim, 3,3, 10,10, domain);
    fit.initSolver(param0,param1);
    fit.solve(values);

    // visualize
    ON_NurbsSurface nurbs = fit.getSurface();
    TomGine::tgRenderModel mesh;
    mesh.SetColor(50,50,200);
    Nurbs2Mesh(nurbs,mesh,domain);
    viewer.AddModel3D(mesh);
  }

  viewer.Update();
  viewer.WaitForEvent(TomGine::TMGL_Press, TomGine::TMGL_Escape);
  return 0;
}
