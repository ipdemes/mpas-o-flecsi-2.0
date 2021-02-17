/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "mpasoflecsi/tasks/io.hh"
#include "mpasoflecsi/specialization/definition.hh"

#include "state.hh"
#include "tasks/testcase.hh"
#include "initialize.hh"

using namespace flecsi;

namespace mpas { namespace sw {

namespace task {

static void init_bottom_depth(mesh::accessor<ro, ro> m,
                              field<double>::accessor<wo, na> h_s) {
  for (auto c : m.cells()) {
    h_s(c) = 4500;
  }
}

}

namespace action {

int init_mesh()
{
  flog(info) << "Initializing mesh: " << mesh_filename.value() << std::endl;

  coloring.allocate(mesh_filename.value());
  m.allocate(coloring.get());

  return 0;
}


int read_fields()
{
  flog(info) << "Reading mesh fields" << std::endl;

  using namespace io::h5;

  hid_t file = open_file(mesh_filename.value(), H5F_ACC_RDONLY, H5P_DEFAULT);

  execute<mpas::task::read_mesh_fields, mpi>(m, file, latCell(m), lonCell(m),
                                             latVertex(m), lonVertex(m),
                                             latEdge(m), lonEdge(m),
                                             xCell(m), yCell(m), zCell(m),
                                             areaCell(m), areaTriangle(m), kiteAreasOnVertex(m),
                                             dvEdge(m), dcEdge(m), meshDensity(m), fVertex(m), weightsOnEdge(m));
  close_file(file);

  execute<task::init_bottom_depth>(m, bottomDepth(m));

  return 0;
}


int init_testcase()
{
  auto tc_num = test_case.value();
  if (tc_num == 1) {
    execute<task::setup_case_1>(m, dvEdge(m), latCell(m), lonCell(m),
                                latVertex(m), lonVertex(m), u(m), h(m));
  }
  return 0;
}

}}}
