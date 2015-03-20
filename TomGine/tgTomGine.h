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

#ifndef _TOMGINE_H_
#define _TOMGINE_H_

#include "tgCamera.h"
#include "tgCollission.h"
#include "tgEngine.h"
#include "tgError.h"
#include "tgErrorMetric.h"
#include "tgFont.h"
#include "tgFrameBufferObject.h"
#include "tgFrustum.h"
#include "tgLight.h"
#include "tgMaterial.h"
#include "tgMathlib.h"
#include "tgModel.h"
#include "tgModelLoader.h"
#include "tgPlot2D.h"
#include "tgPose.h"
#include "tgQuaternion.h"
#include "tgRenderModel.h"
#include "tgShapeCreator.h"
#include "tgTexture.h"
#include "tgTextureModel.h"
#include "tgTimer.h"
#include "tgTomGineThread.h"

/** @brief TomGine namespace */
namespace TomGine{

static void cut_file_name(std::string full_file_name, std::string& file_name, std::string& path)
{
  size_t c_slash_idx = full_file_name.find_last_of("/\\");
  size_t c_dot_idx = full_file_name.find_last_of(".");
  if(c_slash_idx == std::string::npos)
    c_slash_idx = 0;
  else
    c_slash_idx++;
  if(c_dot_idx == std::string::npos || c_dot_idx < c_slash_idx)
    c_dot_idx = full_file_name.size();
  file_name = full_file_name.substr(c_slash_idx,c_dot_idx-c_slash_idx);
  path = full_file_name.substr(0, c_slash_idx);
  if(c_slash_idx==0)
    path= "./";
}

}

#endif
