#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#
set(common_HEADERS
  constants.hh
  fields.hh
  types.hh
)

#------------------------------------------------------------------------------#
# Export header list to parent scope.
#------------------------------------------------------------------------------#

set(common_HEADERS
  ${common_HEADERS}
  PARENT_SCOPE
)


