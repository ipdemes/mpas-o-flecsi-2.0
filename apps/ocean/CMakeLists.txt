#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#
include_directories(${PROJECT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#
add_executable(ocean
  run.cc
  initialize.cc
  tasks/testcase.cc
  tasks/solve.cc
  tasks/init_functions.cc
  ocean.cc
  finalize.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#
target_link_libraries(ocean
  mpasoflecsi
  ${MPAS_LIBRARY_DEPENDENCIES}
)

set (MPASFLECSI_DATA_DIR "${PROJECT_SOURCE_DIR}/data/")
add_test(
  NAME mpasflecsi_ocean
  COMMAND $<TARGET_FILE:ocean> ${MPASFLECSI_DATA_DIR}/meshes/Sphere_240QU.nc -s 5 --num_output 0
)
