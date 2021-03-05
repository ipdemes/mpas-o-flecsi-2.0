/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/mesh.hh>

#include "iohelper.hh"

namespace mpas { namespace task {

void read_mesh_fields(mesh::accessor<flecsi::ro, flecsi::ro> m,
                      hid_t file,
                      io::acc<double> latCell,
                      io::acc<double> lonCell,
                      io::acc<double> latVertex,
                      io::acc<double> lonVertex,
                      io::acc<double> latEdge,
                      io::acc<double> lonEdge,
                      io::acc<double> xCell,
                      io::acc<double> yCell,
                      io::acc<double> zCell,
                      io::acc<double> areaCell,
                      io::acc<double> areaTriangle,
                      io::acc<vdtensor<int>> kiteAreasOnVertex,
                      io::acc<double> dvEdge,
                      io::acc<double> dcEdge,
                      io::acc<double> meshDensity,
                      io::acc<double> fVertex,
                      io::acc<me2tensor<double>> weightsOnEdge);

}}
