/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/control.hh>

#include <flecsi/flog.hh>


namespace mpas { namespace ocean { namespace action {

int finalize();

inline control::action<finalize, cp::finalize> finalize_action;

}}}
