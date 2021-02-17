/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <flecsi/flog.hh>

#include <mpasoflecsi/specialization/control.hh>
#include <mpasoflecsi/specialization/mesh.hh>

#include "state.hh"

namespace mpas {
namespace sw {
namespace action {

int init_mesh();
int read_fields();
int init_testcase();

inline control::action<init_mesh, cp::initialize> init_mesh_action;
inline control::action<read_fields, cp::initialize> read_fields_action;
inline auto const init_dep = read_fields_action.add(init_mesh_action);
inline control::action<init_testcase, cp::initialize> init_testcase_action;
inline auto const tc_dep = init_testcase_action.add(read_fields_action);

}}}
