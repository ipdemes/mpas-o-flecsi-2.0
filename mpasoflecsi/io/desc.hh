/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/common/types.hh"

#include "io_hdf5.hh"
#include "helpers.hh"

namespace mpas { namespace io {

namespace desc {

struct num_vert_levels {
  static constexpr char name[] = "nVertLevels";
};
struct num_tracers {
  static constexpr char name[] = "nTracers";
};


template<mesh::index_space ispace, class T, class... Scales>
struct field {
  using field_type = mpas::acc<T, ro, na>;
  static constexpr mesh::index_space index_space = ispace;
  static constexpr std::size_t ndims = 1 + sizeof...(Scales);
  inline static std::array<const char *, ndims> scales{
    {ispace_size_label<ispace>(), Scales::name...}};
};

struct thickness : field<mesh::index_space::cells,
                         vlreal,
                         num_vert_levels> {
  static constexpr char name[] = "h";
};

struct topography : field<mesh::index_space::cells, double> {
  static constexpr char name[] = "h_s";
};

struct vorticity : field<mesh::index_space::vertices,
                         vlreal,
                         num_vert_levels> {
  static constexpr char name[] = "vorticity";
};

struct pv_vertex : field<mesh::index_space::vertices,
                         vlreal,
                         num_vert_levels> {
  static constexpr char name[] = "potential_vorticity";
};

struct height : field<mesh::index_space::cells, double> {
  static constexpr char name[] = "height";
};

struct f : field<mesh::index_space::vertices, double> {
  static constexpr char name[] = "fVertex";
};

struct tracers : field<mesh::index_space::cells,
                       vltracer,
                       num_vert_levels,
                       num_tracers> {
  static constexpr char name[] = "tracers";
};

}

}}
