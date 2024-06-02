set(SOURCES
    medial_axis_demo.cpp
    polygon_painting.cpp
    polygon_set_painting.cpp
    medial_axis_painting.cpp
    pruned_medial_axis_painting.cpp
    grid_painting.cpp
    pruned_grid_painting.cpp
    visibility_graph_painting.cpp
    backbone_painting.cpp
)

add_executable(medial_axis_demo ${SOURCES})

target_link_libraries(
    medial_axis_demo 
    PRIVATE
    ${COMMON_CLA_TARGET}
    core
    stenomap
    flow_map
    renderer
    CGAL::CGAL
    Qt5::Widgets
    glog::glog
)

install(TARGETS medial_axis_demo DESTINATION ${INSTALL_BINARY_DIR})
