/*
   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */

#include "control.hh"
#include "finalize.hh"
#include "initialize.hh"
#include "options.hh"
#include "remap.hh"
#include "state.hh"

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

int
main(int argc, char ** argv) {
  auto status = flecsi::initialize(argc, argv);
  status = coupler::control::check_status(status);

  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }

  flecsi::log::add_output_stream("clog", std::clog, true);

  status = flecsi::start(coupler::control::execute);

  flecsi::finalize();

  return status;
} // main
