/*
   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

#include "flecsi/data.hh"
#include "flecsi/execution.hh"
#include "flecsi/flog.hh"
#include "flecsi/topo/unstructured/coloring_utils.hh"
#include "flecsi/topo/unstructured/interface.hh"
#include "flecsi/util/parmetis.hh"

#include "mpasoflecsi/common/types.hh"
#include "definition.hh"


namespace mpas {

struct mesh : flecsi::topo::specialization<flecsi::topo::unstructured, mesh> {

  enum index_space { vertices, edges, cells };
  using index_spaces = has<vertices, edges, cells>;
  using connectivities = list<from<cells, to<vertices, edges, cells>>,
                              from<edges, to<vertices, edges, cells>>,
                              from<vertices, to<edges, cells>>>;

  enum entity_list { owned, exclusive, shared, ghost, boundary };
  using entity_lists = list<entity<cells, has<owned, exclusive, shared, ghost>>,
                            entity<edges, has<owned, boundary>>>;

  template<auto>
  static constexpr std::size_t privilege_count = 2;

  template<class B>
  struct interface : B {

    auto cells() {
      return B::template entities<index_space::cells>();
    }

    template<typename B::subspace_list L>
    auto cells() {
      return B::template subspace_entities<index_space::cells, L>();
    }

    template<index_space From>
    auto cells(flecsi::topo::id<From> from) {
      return B::template entities<index_space::cells>(from);
    }

    auto vertices() {
      return B::template entities<index_space::vertices>();
    }

    template<index_space From>
    auto vertices(flecsi::topo::id<From> from) {
      return B::template entities<index_space::vertices>(from);
    }

    auto edges() {
      return B::template entities<index_space::edges>();
    }

    template<index_space From>
    auto edges(flecsi::topo::id<From> from) {
      return B::template entities<index_space::edges>(from);
    }

    template<entity_list List>
    auto edges() {
      return B::template special_entities<index_space::edges, List>();
    }

  }; // struct interface


  struct coloring_policy {
    // primary independent closure token
    using primary =
      flecsi::topo::unstructured_impl::primary_independent<index_space::cells,
        2 /* dimension */,
        0 /* through dimension */,
        1 /* depth */>;

    using auxiliary = std::tuple<
      flecsi::topo::unstructured_impl::auxiliary_independent<index_space::vertices,
        0 /* dimension */,
        2 /* primary dimension */>,
      flecsi::topo::unstructured_impl::auxiliary_independent<index_space::edges,
        1 /* dimension */,
        2 /* primary dimension */>>;

    static constexpr size_t auxiliary_colorings =
      std::tuple_size<auxiliary>::value;
    using definition = io::definition<double>;
  }; // struct coloring_policy


  static coloring color(const std::string & fname);
  static void initialize(flecsi::data::topology_slot<mesh> & s,
                         const coloring & c);


  static void init_cnx(flecsi::field<flecsi::util::id, flecsi::data::ragged>::mutator<flecsi::rw, flecsi::na> e2e,
                       flecsi::topo::unstructured_impl::crs const & cnx);
  static void transpose_cnx(flecsi::field<flecsi::util::id, flecsi::data::ragged>::mutator<flecsi::rw, flecsi::na> v2c,
                            flecsi::field<flecsi::util::id, flecsi::data::ragged>::accessor<flecsi::ro, flecsi::na> c2v);

};


// mesh fields
inline const flecsi::field<int>::definition<mesh, mesh::cells>    maxLevelCell;
inline const flecsi::field<int>::definition<mesh, mesh::edges>    nEdgesOnEdge;

inline const flecsi::field<double>::definition<mesh, mesh::edges>    angleEdge;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    areaCell;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> areaTriangle;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    bottomDepth;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    bottomDepthObserved;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    boundaryLayerDepth;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    dcEdge;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    dvEdge;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    fCell;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    fEdge;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> fVertex;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    latCell;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    latEdge;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> latVertex;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    lonCell;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    lonEdge;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> lonVertex;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    meshDensity;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    surfaceStress;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    xCell;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    xEdge;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> xVertex;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    yCell;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    yEdge;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> yVertex;
inline const flecsi::field<double>::definition<mesh, mesh::cells>    zCell;
inline const flecsi::field<double>::definition<mesh, mesh::edges>    zEdge;
inline const flecsi::field<double>::definition<mesh, mesh::vertices> zVertex;

inline const flecsi::field<vltensor<int>>::definition<mesh, mesh::edges>    edgeMask;

inline const flecsi::field<vltensor<double>>::definition<mesh, mesh::cells>    restingThickness;

inline const flecsi::field<vdtensor<double>>::definition<mesh, mesh::vertices> kiteAreasOnVertex;

inline const flecsi::field<me2tensor<double>>::definition<mesh, mesh::edges> weightsOnEdge;

}
