/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/control.hh>

#include <flecsi/flog.hh>

namespace mpas {
namespace action {

int
init() {

  return 0;
} //init_mesh

control::action<init, cp::initialize> init_action;

} // namespace action
} // namespace mpas
