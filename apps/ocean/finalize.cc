/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/testcase.hh"
#include "state.hh"
#include "finalize.hh"

using namespace flecsi;

namespace mpas { namespace ocean { namespace action {

int finalize() {

  flog(info) << "Finalize run" << std::endl;
  auto elapsed = control::state().elapsed();

  if(inputs::test_case.value() == test::case1) {
    execute<task::calc_error_case_1>(m, latCell(m), lonCell(m), xCell(m), yCell(m), tracers[curr](m), elapsed);
  } else if(inputs::test_case.value() == test::case2) {
    execute<task::calc_error_case_2>(m, latCell(m), lonCell(m), xCell(m), yCell(m), zCell(m), tracers[curr](m), elapsed);
  } else if(inputs::test_case.value() == test::case3) {
    execute<task::calc_error_case_3>(m, latCell(m), lonCell(m), tracers[curr](m), control::state().hmix_coef(), elapsed);
  }

  return 0;
}

}}}
