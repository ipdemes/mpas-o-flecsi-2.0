/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/


#include "remap.hh"
#include "state.hh"
#include "tasks/copy_fields.hh"

#include <flecsi/flog.hh>

using namespace flecsi;

int
coupler::action::remap() {
  execute<task::copy_fields>(m, u(m), v(mm));
  return 0;
} // remap
