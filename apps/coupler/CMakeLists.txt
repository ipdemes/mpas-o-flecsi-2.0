#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#
include_directories(${PROJECT_SOURCE_DIR} include)

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#
add_executable(coupler
  coupler.cc
  initialize.cc
  finalize.cc
  remap.cc
  tasks/copy_fields.cc
  tasks/init_fields.cc
  tasks/io.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#
target_link_libraries(coupler
  mpasoflecsi
  ${MPAS_LIBRARY_DEPENDENCIES}
)
