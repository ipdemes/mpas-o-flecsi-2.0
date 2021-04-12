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
