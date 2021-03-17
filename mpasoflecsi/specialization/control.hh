/*
   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

#include <memory>

#include <flecsi/flog.hh>
#include <flecsi/run/control.hh>

#include <mpasoflecsi/io/state.hh>

namespace mpas {

enum class cp { initialize, run, finalize };

inline const char * operator*(cp control_point) {
  switch(control_point) {
    case cp::initialize:
      return "initialize";
    case cp::run:
      return "run";
    case cp::finalize:
      return "finalize";
  }
  flog_fatal("invalied control point");
}

struct control_policy {

  using control_points_enum = cp;
  struct node_policy {};

  using control = flecsi::run::control<control_policy>;

  template<auto CP>
  using control_point = flecsi::run::control_point<CP>;

  using control_points = std::tuple<control_point<cp::initialize>,
    control_point<cp::run>,
    control_point<cp::finalize>>;

  inline void init_steps(std::size_t num_steps) {
    nsteps = num_steps;
    output_flag = std::make_unique<io::output_flagger>(nsteps);
  }

  inline std::size_t steps() const { return nsteps; }
  inline bool output() const { return (*output_flag); }
  inline bool output_step(std::size_t i) const
  {
    return output_flag->output_step(i);
  }
  inline double & elapsed() { return elapsed_time; }
  inline double elapsed() const { return elapsed_time; }

protected:
  std::size_t nsteps;
  std::unique_ptr<io::output_flagger> output_flag;
  double elapsed_time;
}; // struct control_policy

using control = flecsi::run::control<control_policy>;

} // namespace mpas
