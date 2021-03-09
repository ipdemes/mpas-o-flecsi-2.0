#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#
set(tasks_HEADERS
  io.hh)

set(tasks_SOURCES
  io.cc)

#------------------------------------------------------------------------------#
# Export header list to parent scope.
#------------------------------------------------------------------------------#

set(tasks_HEADERS
  ${tasks_HEADERS}
  PARENT_SCOPE
)

set(tasks_SOURCES
  ${tasks_SOURCES}
  PARENT_SCOPE
)
