/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "control.hh"

namespace coupler {
namespace action {

int finalize();
inline control::action<finalize, cp::finalize> finalize_action;

} // namespace action
} // namespace coupler
