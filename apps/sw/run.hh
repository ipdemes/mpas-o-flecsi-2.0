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
run() {

 return 0;
} //init_mesh

control::action<run, cp::run> run_action;

} // namespace action
} // namespace mpas
