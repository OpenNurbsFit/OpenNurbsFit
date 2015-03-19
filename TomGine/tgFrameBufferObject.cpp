/*
 * Software License Agreement (GNU General Public License)
 *
 *  Copyright (c) 2011, Thomas Mörwald
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @author thomas.moerwald
 *
 */
#include "tgFrameBufferObject.h"
#include "tgError.h"
#include "tgMathlib.h"
//#include <opencv2/highgui/highgui.hpp>
#include <stdexcept>

using namespace TomGine;

tgFrameBufferObject::tgFrameBufferObject(unsigned w, unsigned h,
                                         GLint colorInternal, GLint depthInternal)
{

  m_width = w;
  m_height = h;

  texColor.assign(1, new tgTexture2D());
  texDepth = new tgTexture2D();

  texColor[0]->Bind();
  texColor[0]->Load(NULL, m_width, m_height, colorInternal, GL_RGBA, GL_UNSIGNED_BYTE);
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] fbo_tex");

  texDepth->Bind();
  texDepth->Load(NULL, m_width, m_height, depthInternal, GL_DEPTH_COMPONENT, GL_FLOAT);
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] fbo_depth_tex");

  glGenFramebuffers(1, &m_fbo_id);
  Bind();

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColor[0]->GetTextureID(), 0);
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] attach color texture");

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, texDepth->GetTextureID(), 0);
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] attach depth texture");

  Unbind();

  if (tgCheckFBError(GL_FRAMEBUFFER, "[tgFrameBufferObject::tgFrameBufferObject]") != GL_FRAMEBUFFER_COMPLETE
      || tgCheckError("[tgFrameBufferObject::tgFrameBufferObject]") != GL_NO_ERROR)
  {
    std::string errmsg =
        std::string("[tgFrameBufferObject::tgFrameBufferObject] Error generating frame buffer objects");
    throw std::runtime_error(errmsg.c_str());
  }
  glDisable(GL_TEXTURE_2D);
}

tgFrameBufferObject::tgFrameBufferObject(unsigned w, unsigned h,
                                         std::vector<GLint>& colorInternal, GLint depthInternal)
{

  m_width = w;
  m_height = h;

  for(size_t i=0; i<colorInternal.size(); i++)
  {
    texColor.push_back(new tgTexture2D());
    texColor[i]->Bind();
    texColor[i]->Load(NULL, m_width, m_height, colorInternal[i], GL_RGBA, GL_UNSIGNED_BYTE);
  }
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] fbo_tex");

  texDepth = new tgTexture2D();
  texDepth->Bind();
  texDepth->Load(NULL, m_width, m_height, depthInternal, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE);
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] fbo_depth_tex");

  glGenFramebuffers(1, &m_fbo_id);
  Bind();

  std::vector<GLenum> drawBuffers;
  for(size_t i=0; i<texColor.size(); i++)
  {
    drawBuffers.push_back(GL_COLOR_ATTACHMENT0+i);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0+i, GL_TEXTURE_2D, texColor[i]->GetTextureID(), 0);
  }
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] attach color texture");

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, texDepth->GetTextureID(), 0);
  tgCheckError("[tgFrameBufferObject::tgFrameBufferObject] attach depth texture");

  glDrawBuffers(drawBuffers.size(), &drawBuffers[0]);

  Unbind();

  if (tgCheckFBError(GL_FRAMEBUFFER, "[tgFrameBufferObject::tgFrameBufferObject]") != GL_FRAMEBUFFER_COMPLETE
      || tgCheckError("[tgFrameBufferObject::tgFrameBufferObject]") != GL_NO_ERROR)
  {
    std::string errmsg =
        std::string("[tgFrameBufferObject::tgFrameBufferObject] Error generating frame buffer objects");
    throw std::runtime_error(errmsg.c_str());
  }
  glDisable(GL_TEXTURE_2D);
}

tgFrameBufferObject::~tgFrameBufferObject()
{
  if (glIsFramebuffer(m_fbo_id))
    glDeleteFramebuffers(1, &m_fbo_id);

  for(size_t i=0; i<texColor.size(); i++)
    delete texColor[i];
  delete texDepth;
}

void tgFrameBufferObject::Bind()
{
  glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_id);
  tgCheckError("[tgFrameBufferObject::Bind]");
}

void tgFrameBufferObject::Unbind()
{
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void tgFrameBufferObject::GetPointCloud(cv::Mat4f &points, bool world_coords)
{
  // get values of rendering pipeline
  mat4 proj, model;
  GLint view[4];
  glGetFloatv(GL_PROJECTION_MATRIX, proj);
  glGetFloatv(GL_MODELVIEW_MATRIX, model);
  glGetIntegerv(GL_VIEWPORT, view);

  cv::Vec4f cvp;
  cvp[0] = NAN;
  cvp[1] = NAN;
  cvp[2] = NAN;
  cvp[3] = NAN;
  points.setTo(cvp);

  float depth_range[2];
  glGetFloatv(GL_DEPTH_RANGE, depth_range);

  mat3 Rt(model.transpose());
  vec3 t(model[12], model[13], model[14]);

  const float &fx = proj[0];
  const float &fy = proj[5];
  const float &cx = proj[8];
  const float &cy = proj[9];
  const float &z1 = proj[10];
  const float &z2 = proj[14];
  const float &dNear = depth_range[0];
  const float &dFar = depth_range[1];
  const int &width = view[2];
  const int &height = view[3];

  float zFar = z2 / (z1 + 1.0f);
  float zNear = z2 * zFar / (z2 - 2.0f * zFar);

  cv::Mat4b color(m_height,m_width);
  cv::Mat1f depth(m_height,m_width);

  texColor[0]->Bind();
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, color.data);

  texDepth->Bind();
  glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_FLOAT, depth.data);

  // parse depth values and create 3D points in world coordinates
  for (int j = 0; j < depth.rows; j++)
  {
    for (int i = 0; i < depth.cols; i++)
    {
      int winY = j;
      if (world_coords)
        winY = height - j - 1; // inverse y coordinates

      int idx = (depth.rows - j - 1) * depth.cols + i; // depthbuffer index
      const float &z_b = depth(idx);
      const cv::Vec4b &c = color(idx);
      RGBValue col;
      col.r = c(0);
      col.g = c(1);
      col.b = c(2);
      col.a = c(3);


      if (z_b > dNear && z_b < dFar) // depth range check
      {

        float z_n = 2.0 * z_b - 1.0; // normalized depth

        // transform to camera coordinates
        float z = 2.0 * zNear * zFar / (zFar + zNear - z_n * (zFar - zNear));
        float x = (2.0 * float(i) / width - 1.0f + cx) * z / fx;
        float y = (2.0 * float(winY) / height - 1.0f + cy) * z / fy;

        vec4 p(x, y, z, 1);
        if (world_coords) // transform to world coordinates
        {
          p.z = -p.z;
          p = Rt * (p - t);
        }

        cvp[0] = p.x;
        cvp[1] = p.y;
        cvp[2] = p.z;
        cvp[3] = col.float_value;

        points(j, i) = cvp;

      }
    }
  }
}

//void tgFrameBufferObject::SaveColor(const char* filename)
//{
//  texColor[0]->Bind();
//  cv::Mat img(m_height, m_width, CV_8UC3);
//  glGetTexImage(GL_TEXTURE_2D, 0, GL_BGR, GL_UNSIGNED_BYTE, img.data);
//  glDisable(GL_TEXTURE_2D);
//  tgCheckError("[tgFrameBufferObject::SaveColor]");
//  cv::imwrite(filename, img);
//}

//void tgFrameBufferObject::SaveDepth(const char* filename)
//{
//  texDepth->Bind();
//  cv::Mat img(m_height, m_width, CV_8U);
//  glGetTexImage(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, img.data);
//  glDisable(GL_TEXTURE_2D);
//  tgCheckError("[tgFrameBufferObject::SaveDepth]");
//  cv::imwrite(filename, img);
//}
