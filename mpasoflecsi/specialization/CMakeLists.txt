#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#

set(specialization_HEADERS
  control.hh
  mesh.hh)

set(specialization_SOURCES
  mesh.cc)

#------------------------------------------------------------------------------#
# Export header list to parent scope.
#------------------------------------------------------------------------------#

set(specialization_HEADERS
  ${specalization_HEADERS}
  PARENT_SCOPE
)

set(specialization_SOURCES
  ${specialization_SOURCES}
  PARENT_SCOPE)
