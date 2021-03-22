/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "finalize.hh"
#include "control.hh"
#include "state.hh"
#include "tasks/io.hh"

#include <flecsi/flog.hh>

using namespace flecsi;


int
coupler::action::finalize() {
  execute<task::io>(mesh_src, mesh_trg, ud(mesh_src), vd(mesh_trg));
  return 0;
} // remap

