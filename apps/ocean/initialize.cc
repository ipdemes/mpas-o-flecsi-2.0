/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "mpasoflecsi/tasks/io.hh"
#include "mpasoflecsi/specialization/definition.hh"

#include "state.hh"
#include "tasks/testcase.hh"
#include "tasks/init_functions.hh"
#include "initialize.hh"

using namespace flecsi;

namespace mpas { namespace ocean { namespace action {


int init_mesh()
{
  flog(info) << "Initializing mesh: " << inputs::meshfile.value() << std::endl;

  coloring.allocate(inputs::meshfile.value());
  m.allocate(coloring.get());

  return 0;
}


int read_fields()
{
  flog(info) << "Reading mesh fields" << std::endl;

  using namespace io::h5;

  hid_t file = open_file(inputs::meshfile.value(), H5F_ACC_RDONLY, H5P_DEFAULT);

  execute<mpas::task::read_ocean_mesh_fields, mpi>(m, file, angleEdge(m), areaCell(m), areaTriangle(m),
                                             bottomDepth(m), bottomDepthObserved(m), dcEdge(m),
                                             dvEdge(m), edgeMask(m), fCell(m), fEdge(m), fVertex(m),
                                             kiteAreasOnVertex(m), latCell(m), latEdge(m), latVertex(m),
                                             //layerThickness(m), 
                                             lonCell(m), lonEdge(m), lonVertex(m), //maxLevelCell(m),
                                             meshDensity(m), nEdgesOnEdge(m),// normalVelocity[curr](m),
                                             restingThickness(m), tracers[curr](m), weightsOnEdge(m),
                                             xCell(m), xEdge(m), xVertex(m), yCell(m), yEdge(m),
                                             yVertex(m), zCell(m), zEdge(m), zVertex(m)
                                             );

  close_file(file);

  execute<task::init_extra_fields>(m, w(m), boundaryCell(m), boundaryLayerDepth(m), 
                                  surfaceStress(m), maxLevelCell(m));

  return 0;
}


int init_testcase()
{
  control::state().init_steps(inputs::nsteps.value());
  flog(info) << "Initialize test case" << std::endl;

  switch (inputs::test_case.value()) {
  case test::case1 :
    control::state().rk4_config() = true;
    control::state().monoAdvect_config() = true;
    execute<task::setup_case_1>(m, angleEdge(m), latEdge(m), xCell(m), normalVelocity[curr](m),
                                layerThickness[curr](m), tracers[curr](m));
    break;
  case test::case2 :
    control::state().rk4_config() = true;
    control::state().monoAdvect_config() = true;
    execute<task::setup_case_2>(m, angleEdge(m), latEdge(m), latCell(m), lonCell(m), xCell(m), yCell(m),
                                zCell(m), normalVelocity[curr](m), layerThickness[curr](m), tracers[curr](m));
  case test::case3 :
    control::state().rk4_config() = true;
    control::state().hmix_config() = true;
    control::state().hmix_coef() = 1e7;
    execute<task::setup_case_3>(m, angleEdge(m), latEdge(m), latCell(m), lonCell(m), normalVelocity[curr](m),
                                layerThickness[curr](m), tracers[curr](m), control::state().hmix_coef());
    break;
  }

  return 0;
}

int derive_initial_fields()
{
  flog(info) << "Derive remaining fields" << std::endl;

  execute<task::init_maxLevelEdgeTop>(m, maxLevelCell(m), maxLevelEdgeTop(m));
  execute<task::init_layerThicknessEdge>(m, layerThickness[curr](m), maxLevelEdgeTop(m),
                                         layerThicknessEdge(m));
  execute<task::init_meshScalingDel2>(m, meshDensity(m), meshScalingDel2(m));
  execute<task::init_edgeSignOnCell>(m, edgeSignOnCell(m));
  execute<task::init_2nd_deriv>(m, angleEdge(m), dcEdge(m), xCell(m), yCell(m), zCell(m),
                                xVertex(m), yVertex(m), zVertex(m), derivTwo(m));
  execute<task::init_tracer_adv_coeff>(m, maxLevelCell(m), dcEdge(m), dvEdge(m), boundaryCell(m),
                                       derivTwo(m), nAdvCellsForEdge(m), highOrderAdvMask(m),
                                       advCellsForEdge(m), advCoef(m), advCoef_3rd(m));

  return 0;
}

}}}
