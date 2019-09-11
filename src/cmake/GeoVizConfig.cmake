include(CMakeFindDependencyMacro)

# Capturing values from configure (optional)
#set(my-config-var @my-config-var@)

# Same syntax as find_package
#find_dependency(MYDEP REQUIRED)

# Any extra setup

# Add the targets file
include("${CMAKE_CURRENT_LIST_DIR}/GeoVizTargets.cmake")