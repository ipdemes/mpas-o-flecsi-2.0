/*
 *    Copyright (c) 2016, Triad National Security, LLC
 *       All rights reserved.
 *                                                                                     */

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include <mpasoflecsi/specialization/control.hh>

#include "finalize.hh"
#include "initialize.hh"
#include "run.hh"

int main(int argc, char *argv[])
{
  auto status = flecsi::initialize(argc, argv);

  status = mpas::control::check_status(status);

  if(status != flecsi::run::status::success) {
    return status < flecsi::run::status::clean ? 0 : status;
  }

  flecsi::log::add_output_stream("clog", std::clog, true);

  status = flecsi::start(mpas::control::execute);

  flecsi::finalize();

  return status;
}
