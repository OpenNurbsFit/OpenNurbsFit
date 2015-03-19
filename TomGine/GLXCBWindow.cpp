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

// http://renderingpipeline.com/2012/05/windowless-opengl/

#include "GLWindow.h"
#include <stdio.h>
#include <stdexcept>
#include <vector>

namespace TomGine {

GLWindow::GLWindow() : windowless(true)
{
  //http://sidvind.com/wiki/Opengl/windowless
  GLint attr[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };

  // open display
  if ( ! (dpy = XOpenDisplay(NULL)) ) {
    fprintf(stderr, "cannot connect to X server\n\n");
      exit(1);
  }
  // get root window
  root = DefaultRootWindow(dpy);

  // get visual matching attr
  if( ! (vi = glXChooseVisual(dpy, 0, attr)) ) {
       fprintf(stderr, "no appropriate visual found\n\n");
       exit(1);
  }
  // create a context using the root window
  if ( ! (glc = glXCreateContext(dpy, vi, NULL, GL_TRUE)) ){
      fprintf(stderr, "failed to create context\n\n");
      exit(1);
  }
  glXMakeCurrent(dpy, root, glc);
  // try it out, remember to *NOT* render to the default framebuffer!
}

GLWindow::GLWindow(unsigned int width, unsigned int height, const char* name, bool threaded,
                   bool stereo) : windowless(false)
{
  XInitThreads();

  dpy = XOpenDisplay(NULL);

  if (dpy == NULL) {
    throw std::runtime_error("[GLWindow::init] Error cannot connect to X server");
  }

  root = DefaultRootWindow(dpy);

  if (stereo) {
    GLint att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, GLX_STEREO, None };
    vi = glXChooseVisual(dpy, 0, att);
  } else {
    GLint att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
    vi = glXChooseVisual(dpy, 0, att);
  }

  if (vi == NULL)
    throw std::runtime_error("[GLWindow::init] Error no appropriate visual found");

  cmap = XCreateColormap(dpy, root, vi->visual, AllocNone);

  swa.colormap = cmap;
  swa.event_mask = ExposureMask | KeyPressMask | KeyReleaseMask | ButtonPressMask
      | ButtonReleaseMask | PointerMotionMask | StructureNotifyMask;

  glWin = XCreateWindow(dpy, root, 0, 0, width, height, 0, vi->depth, InputOutput, vi->visual,
                        CWColormap | CWEventMask, &swa);
  wmDelete = XInternAtom(dpy, "WM_DELETE_WINDOW", true);
  XSetWMProtocols(dpy, glWin, &wmDelete, 1);

  XMapWindow(dpy, glWin);
  XStoreName(dpy, glWin, name);

  glc = glXCreateContext(dpy, vi, NULL, GL_TRUE);
  glXMakeCurrent(dpy, glWin, glc);

  connection = XGetXCBConnection(dpy);
  if(!connection) {
    XCloseDisplay(dpy);
    fprintf(stderr, "Can't get xcb connection from display\n");
    exit(-1);
  }

  printf("GLWindow '%s' %dx%d\nOpenGL Version: %s\n", name, width, height, glGetString(GL_VERSION));

  /* Load the font. */
#ifdef GLX_FONT
  font_info = XLoadQueryFont(dpy, GLX_FONT);
  font_base = glGenLists(256);
  if (!font_info) {
    fprintf(stderr, "XLoadQueryFont() failed - Exiting.\n");
    exit(-1);
  }
  else {
    /* Tell GLX which font & glyphs to use. */
    int first = font_info->min_char_or_byte2;
    int last = font_info->max_char_or_byte2;
    glXUseXFont(font_info->fid, first, last-first+1, font_base+first);
  }
#endif

  this->threaded = threaded;

  glXSwapBuffers(dpy, glWin);
}

GLWindow::~GLWindow()
{
  glXMakeCurrent(dpy, None, NULL);
#ifdef GLX_FONT
  if(!windowless)
  {
    glDeleteLists(font_base, 256);
    XFreeFont(dpy, font_info);
  }
#endif
  glXDestroyContext(dpy, glc);
  if(!windowless)
  {
    XFlush(dpy);
    XDestroyWindow(dpy, glWin);
    XFreeColormap(dpy, cmap);
  }
  XCloseDisplay(dpy);
}

void GLWindow::Activate()
{
  if(windowless)
    glXMakeCurrent(dpy, root, glc);
  else
    glXMakeCurrent(dpy, glWin, glc);
}

void GLWindow::Update()
{
  if(!windowless)
    glXSwapBuffers(dpy, glWin);
}

} /* namespace */
