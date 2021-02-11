/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <mpasoflecsi/specialization/mesh.hh>

#include "iohelper.hh"

namespace mpas { namespace task {

template<class T>
using io_acc = io::acc<T>;

void read_mesh_fields(mesh::accessor<flecsi::ro, flecsi::ro> m,
                      hid_t file,
                      io_acc<double> latCell,
                      io_acc<double> lonCell,
                      io_acc<double> latVertex,
                      io_acc<double> lonVertex,
                      io_acc<double> latEdge,
                      io_acc<double> lonEdge,
                      io_acc<double> xCell,
                      io_acc<double> yCell,
                      io_acc<double> zCell,
                      io_acc<double> areaCell,
                      io_acc<double> areaTriangle,
                      io_acc<vdtensor<int>> kiteAreasOnVertex,
                      io_acc<double> dvEdge,
                      io_acc<double> dcEdge,
                      io_acc<double> meshDensity,
                      io_acc<double> fVertex,
                      io_acc<metensor<double>> weightsOnEdge);

}}
