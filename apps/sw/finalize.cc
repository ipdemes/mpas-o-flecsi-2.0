/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/solve.hh"
#include "state.hh"
#include "finalize.hh"

using namespace flecsi;

namespace mpas { namespace sw {
namespace action {

int finalize()
{
  if (inputs::test_case.value() == test::case1) {
    flog(info) << "Computing L2 error" << std::endl;

    auto elapsed = control::state().elapsed();
    auto err_fut = execute<task::compute_error>(m, h[curr](m), latCell(m), lonCell(m), areaCell(m),
                                                elapsed);
    auto err = err_fut.get();
    double tol = 0.0015;
    if (err > tol) throw;
  }

  return 0;
}

}
}}
