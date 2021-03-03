/*
   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

#include <flecsi/flog.hh>
#include <flecsi/run/control.hh>

namespace coupler {

enum class cp { initialize, remap, finalize };

inline const char * operator*(cp control_point) {
  switch(control_point) {
    case cp::initialize:
      return "initialize";
    case cp::remap:
      return "remap";
    case cp::finalize:
      return "finalize";
  }
  flog_fatal("invalid control point");
}

struct control_policy {

  using control_points_enum = cp;
  struct node_policy {};

  using control = flecsi::run::control<control_policy>;

  template<auto CP>
  using control_point = flecsi::run::control_point<CP>;

  using control_points = std::tuple<control_point<cp::initialize>,
    control_point<cp::remap>,
    control_point<cp::finalize>>;
}; // struct control_policy

using control = flecsi::run::control<control_policy>;

} // namespace coupler
