/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "io.hh"

namespace mpas { namespace task {

using namespace flecsi;

void read_mesh_fields(mesh::accessor<ro, ro> m,
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
                      io_acc<metensor<double>> weightsOnEdge)
{
  using namespace mpas::io;

  constexpr auto cell_space = mesh::index_space::cells;
  constexpr auto edge_space = mesh::index_space::edges;
  constexpr auto vert_space = mesh::index_space::vertices;
  { // read 1D float datasets
    helper1d<double> reader(file);
    reader.read<vert_space>("fVertex", fVertex, m);
    reader.read<cell_space>("areaCell", areaCell, m);
    reader.read<vert_space>("areaTriangle", areaTriangle, m);
    reader.read<edge_space>("dvEdge", dvEdge, m);
    reader.read<edge_space>("dcEdge", dcEdge, m);
    reader.read<cell_space>("meshDensity", meshDensity, m);
    reader.read<cell_space>("xCell", xCell, m);
    reader.read<cell_space>("yCell", yCell, m);
    reader.read<cell_space>("zCell", zCell, m);
    reader.read<cell_space>("latCell", latCell, m);
    reader.read<cell_space>("lonCell", lonCell, m);
    reader.read<vert_space>("latVertex", latVertex, m);
    reader.read<vert_space>("lonVertex", lonVertex, m);
    reader.read<vert_space>("latEdge", latEdge, m);
    reader.read<vert_space>("lonEdge", lonEdge, m);
  }
}

}}
