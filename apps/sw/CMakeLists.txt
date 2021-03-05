#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#
include_directories(${PROJECT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#
add_executable(sw
  run.cc
  initialize.cc
  tasks/testcase.cc
  tasks/solve.cc
  sw.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#
target_link_libraries(sw
  mpasoflecsi
  ${MPAS_LIBRARY_DEPENDENCIES}
)

set (MPASFLECSI_DATA_DIR "${PROJECT_SOURCE_DIR}/data/")
add_test(
  NAME mpasflecsi_sw
  COMMAND $<TARGET_FILE:sw> ${MPASFLECSI_DATA_DIR}/meshes/Sphere_240QU.nc
)
