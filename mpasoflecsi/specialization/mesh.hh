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

#include "definition.hh"


namespace mpas {

struct mesh : flecsi::topo::specialization<flecsi::topo::unstructured, mesh> {

  enum index_space { vertices, edges, cells };
  using index_spaces = has<vertices, edges, cells>;
  using connectivities = list<from<cells, to<vertices, edges, cells>>,
                              from<edges, to<vertices, edges, cells>>,
                              from<vertices, to<edges, cells>>>;

  enum entity_list { boundary };
  using entity_lists = list<entity<edges, has<boundary>>>;

  template<auto>
  static constexpr std::size_t privilege_count = 2;

  template<class B>
  struct interface : B {

    auto cells() {
      return B::template entities<index_space::cells>();
    }

    template<index_space From>
    auto cells(flecsi::topo::id<From> from) {
      return B::template entities<index_space::cells>(from);
    }

    auto vertices() {
      return B::template entities<index_space::vertices>();
    }

    auto edges() {
      return B::template entities<index_space::edges>();
    }

    template<entity_list List>
    auto edges() {
      return B::template special_entities<mesh::edges, List>();
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

}; // struct mpas_mesh

} // namespace poisson
