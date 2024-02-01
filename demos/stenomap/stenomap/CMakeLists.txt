set(SOURCES
    stenomap_demo.cpp
)

add_executable(stenomap_demo ${SOURCES})

target_link_libraries(
    stenomap_demo
    PRIVATE
    ${COMMON_CLA_TARGET}
    core
    flow_map
    renderer
    CGAL::CGAL
    Qt5::Widgets
    glog::glog
)

install(TARGETS stenomap_demo DESTINATION ${INSTALL_BINARY_DIR})
