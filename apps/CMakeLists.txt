# find dependent packages
if(WIN32 AND NOT CYGWIN)
  set(EIGEN3_INCLUDE_DIRS "C:/Eigen/" CACHE PATH "Include directory for Eigen3 library.")
else()
  find_package(PkgConfig REQUIRED)
  pkg_search_module(EIGEN3 REQUIRED eigen3)
  include_directories(${EIGEN3_INCLUDE_DIRS})
endif()

if(SUITESPARSE_FOUND)
  include_directories(${SUITESPARSE_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Error, suitesparse not found.")
endif()

# example_surface_depth
add_executable(fit_patch example_fitpatch.cpp triangulation.h triangulation.cpp)
target_link_libraries(fit_patch nurbsfit opennurbs tomgine)

add_executable(fit_opencurve2d example_opencurve2d.cpp triangulation.h triangulation.cpp)
target_link_libraries(fit_opencurve2d nurbsfit opennurbs tomgine)
 
install(TARGETS fit_patch fit_opencurve2d
        EXPORT OpenNurbsFitTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
