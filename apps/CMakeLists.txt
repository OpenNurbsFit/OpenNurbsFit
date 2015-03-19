
# example_surface_depth
add_executable(surface_depth example_surface_depth.cpp)
target_link_libraries(surface_depth nurbsfit opennurbs tomgine)
 
install(TARGETS surface_depth
        EXPORT OpenNurbsFitTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
