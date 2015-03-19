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

# source files
set(SOURCE_CPP
    surface.cpp
    patch.cpp)

set(SOURCE_H
    surface.h
    patch.h)

# build
add_library(nurbsfit SHARED ${SOURCE_CPP} ${SOURCE_H})
target_link_libraries(nurbsfit opennurbs ${SUITESPARSE_LIBRARIES})

# install
install(FILES ${SOURCE_H} DESTINATION "${INSTALL_INCLUDE_DIR}/NurbsFit")
install(TARGETS nurbsfit
        EXPORT OpenNurbsFitTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/NurbsFit"
        COMPONENT dev)

