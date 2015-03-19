set(V4R_TOMGINE 1 CACHE BOOL "OpenGL-GPU rendering, geometry and image processing library." )
IF(V4R_TOMGINE)

set(TG_NO_EVENT_HANDLING 0 CACHE BOOL "Disable threaded event handling (workaround GUI update bug Ubuntu 14.04)")
IF(TG_NO_EVENT_HANDLING)
add_definitions(-DTG_NO_EVENT_HANDLING)
endif(TG_NO_EVENT_HANDLING)

add_definitions(-DLINUX)
add_definitions(-DGL_GLEXT_PROTOTYPES)

# find OpenGL
find_package(PkgConfig REQUIRED)
pkg_search_module(GL REQUIRED gl)
pkg_search_module(GLU REQUIRED glu)

# find OpenCV
find_package( OpenCV REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})

if(NOT GL_FOUND)
	message(FATAL "libgl not found, required by v4rTomGine")
endif(NOT GL_FOUND)

# font selection
set(ENABLE_FTGL_FONTS 0 CACHE BOOL "Enable ftgl fonts do draw OpenGL fonts?" )
IF(ENABLE_FTGL_FONTS)
	pkg_search_module(FTGL REQUIRED ftgl)
	if(NOT FTGL_FOUND)
		message(FATAL "libftgl not found, required by v4rTomGine")
	endif(NOT FTGL_FOUND)
  set(TRUETYPE_FONT /usr/share/fonts/truetype/freefont/FreeSans.ttf CACHE PATH "A valid truetype font on your system")
  add_definitions(-DTTF_FONT="${TRUETYPE_FONT}")
  include_directories(${FTGL_INCLUDE_DIRS})
  add_definitions(-DUSE_FTGL_FONT)
ELSE(ENABLE_FTGL_FONTS)
  set(GLX_FONT_NAME 10x20 CACHE PATH "A valid system font (shell command: xlsfonts)")
  add_definitions(-DGLX_FONT="${GLX_FONT_NAME}")
ENDIF(ENABLE_FTGL_FONTS)

# create tomgine include directory name
string(REGEX REPLACE "v4r" "" TOMGINE_INCLUDE_NAME tomgine)

set(SOURCE_CPP 
  GLXCBWindow.cpp
  GLXCBWindowGetEvent.cpp
	ply.c
	tgCamera.cpp
  tgCollission.cpp
	tgEngine.cpp
	tgError.cpp
	tgErrorMetric.cpp
	tgFont.cpp
	tgFrameBufferObject.cpp
	tgFrustum.cpp
	tgLight.cpp
	tgMaterial.cpp
	tgMathlib.cpp
	tgModel.cpp
	tgModelLoader.cpp
	tgPlot2D.cpp
	tgPose.cpp
	tgQuaternion.cpp
	tgRenderModel.cpp
	tgShapeCreator.cpp
	tgTexture.cpp
	tgTextureModel.cpp
	tgTimer.cpp
	tgTomGineThread.cpp)

set(SOURCE_H
	GLEvent.h
	GLInput.h
	GLWindow.h
	headers.h
	ply.h
	PlyStructure.h
	tgCamera.h
  tgCollission.h
	tgEngine.h
	tgError.h
	tgErrorMetric.h
	tgFont.h
	tgFrameBufferObject.h
	tgFrustum.h
	tgLight.h
	tgMaterial.h
	tgMathlib.h
	tgModel.h
	tgModelLoader.h
	tgPlot2D.h
	tgPose.h
	tgQuaternion.h
	tgRenderModel.h
	tgShapeCreator.h
	tgSingleton.h
	tgTexture.h
	tgTextureModel.h
	tgTimer.h
	tgTomGine.h
	tgTomGineThread.h
	ub3Image.h
)

# build
add_library(tomgine SHARED ${SOURCE_CPP} ${SOURCE_H})
target_link_libraries(tomgine ${OpenCV_LIBS} ${GL_LIBRARIES} ${GLU_LIBRARIES} pthread X11 X11-xcb xcb rt)
if(ENABLE_FTGL_FONTS)
target_link_libraries(tomgine ${FTGL_LIBRARIES})
endif(ENABLE_FTGL_FONTS)

# install
install(FILES ${SOURCE_H} DESTINATION "${INSTALL_INCLUDE_DIR}/TomGine")
install(TARGETS tomgine
        EXPORT OpenNurbsFitTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/TomGine"
        COMPONENT dev)

ENDIF(V4R_TOMGINE)
