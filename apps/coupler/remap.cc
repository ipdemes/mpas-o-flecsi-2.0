/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/


#include "remap.hh"
#include "state.hh"
#include "tasks/copy_fields.hh"
#include "tasks/init_fields.hh"

#include <flecsi/execution.hh>

using namespace flecsi;

int
coupler::action::remap() {
  execute<task::init_fields>(m, ud(m), vd(m));
  execute<task::copy_fields>(m, ud(m), vd(m));
  return 0;
} // remap
