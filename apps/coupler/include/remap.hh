/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "control.hh"

namespace coupler {
namespace action {

int remap();
inline control::action<remap, cp::remap> remap_action;

} // namespace action
} // namespace coupler
