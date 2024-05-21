set(SOURCES
	medial_axis.cpp
	backbone.cpp
)
set(HEADERS
	medial_axis.h
	backbone.h
)

add_library(stenomap ${SOURCES})
target_link_libraries(stenomap
	PRIVATE glog::glog
)

cartocrow_install_module(stenomap)
install(FILES ${HEADERS} DESTINATION ${CARTOCROW_INSTALL_DIR}/stenomap)
