#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#
include_directories(${PROJECT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#
add_executable(sw
  initialize.cc
  sw.cc
  tasks/testcase.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#
target_link_libraries(sw
  mpasoflecsi
  ${MPAS_LIBRARY_DEPENDENCIES}
)
