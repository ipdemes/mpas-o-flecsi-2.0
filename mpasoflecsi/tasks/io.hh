/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/specialization/mesh.hh"
#include "mpasoflecsi/io/helpers.hh"

namespace mpas { namespace task {

namespace detail {
template<class D, class T>
void
output_field(mesh::accessor<ro, ro> m,
             hid_t file,
             int time_ind,
             T field) {
  typename io::helper<D::ndims, double>::type writer(file);
  writer.template write_timestep<D::index_space>(
    D::name, field, m, time_ind);

  using namespace mpas::io::h5;
  if(time_ind == 0) {
    attach_scale(file, D::name, "Time", 0);
    for(std::size_t i{0}; i < D::scales.size(); ++i) {
      attach_scale(file, D::name, D::scales[i], i + 1);
    }
  }
}
}


template<class D>
struct field_output_wrapper;
template<class... D>
struct field_output_wrapper<type_list<D...>> {
  static void output(mesh::accessor<ro, ro> m,
                     hid_t file,
                     int time_ind,
                     typename D::field_type... fields) {
    (detail::output_field<D>(m, file, time_ind, fields), ...);
  }
};


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
                      io::acc<vdtensor<double>> kiteAreasOnVertex,
                      io::acc<double> dvEdge,
                      io::acc<double> dcEdge,
                      io::acc<double> meshDensity,
                      io::acc<double> fVertex,
                      io::acc<me2tensor<double>> weightsOnEdge);

void read_ocean_mesh_fields(mesh::accessor<flecsi::ro, flecsi::ro> m,
                      hid_t file,
                      io::acc<double> angleEdge,
                      io::acc<double> areaCell,
                      io::acc<double> areaTriangle,
                      io::acc<double> bottomDepth,
                      io::acc<double> bottomDepthObserved,
                      io::acc<double> dcEdge,
                      io::acc<double> dvEdge,
                      io::acc<vltensor<int>> edgeMask,
                      io::acc<double> fCell,
                      io::acc<double> fEdge,
                      io::acc<double> fVertex,
                      io::acc<vdtensor<double>> kiteAreasOnVertex,
                      io::acc<double> latCell,
                      io::acc<double> latEdge,
                      io::acc<double> latVertex,
//                      io::acc<vltensor<double>> layerThickness,
                      io::acc<double> lonCell,
                      io::acc<double> lonEdge,
                      io::acc<double> lonVertex,
//                      io::acc<int> maxLevelCell,
                      io::acc<double> meshDensity,
                      io::acc<int> nEdgesOnEdge,
//                      io::acc<vltensor<double>> normalVelocity,
                      io::acc<vltensor<double>> restingThickness,
                      io::acc<vltracer> tracers,
                      io::acc<me2tensor<double>> weightsOnEdge,
                      io::acc<double> xCell,
                      io::acc<double> xEdge,
                      io::acc<double> xVertex,
                      io::acc<double> yCell,
                      io::acc<double> yEdge,
                      io::acc<double> yVertex,
                      io::acc<double> zCell,
                      io::acc<double> zEdge,
                      io::acc<double> zVertex
                      );


void mesh_output_init(mesh::accessor<ro, ro> m,
                      acc<double, ro, na> areaCell,
                      acc<double, ro, na> dvEdge,
                      acc<double, ro, na> dcEdge,
                      acc<double, ro, na> xCell,
                      acc<double, ro, na> yCell,
                      acc<double, ro, na> zCell,
                      acc<double, ro, na> latCell,
                      acc<double, ro, na> lonCell,
                      acc<double, ro, na> latVertex,
                      acc<double, ro, na> lonVertex,
                      acc<double, ro, na> meshDensity,
                      const char * mfile_name,
                      hid_t * mfile);


void mesh_output_finalize(hid_t mfile);

}}
