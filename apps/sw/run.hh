/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/control.hh>

#include <flecsi/flog.hh>

namespace mpas { namespace sw { namespace action {

int run();

inline control::action<run, cp::run> run_action;

}}}

