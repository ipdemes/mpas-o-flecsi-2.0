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
  using index_spaces = has<cells, vertices, edges>;
  using connectivities = list<from<cells, to<vertices, edges>>>;

  enum entity_list { boundary };
  using entity_lists = list<entity<edges, has<boundary>>>;


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


#if 0
  struct coloring_policy {
    // primary independent closure token
    using primary =
      topo::unstructured_impl::primary_independent<index_space::cells,
        2 /* dimension */,
        0 /* through dimension */,
        1 /* depth */>;

    using auxiliary = std::tuple<
      topo::unstructured_impl::auxiliary_independent<index_space::vertices,
        0 /* dimension */,
        2 /* primary dimension */>,
      topo::unstructured_impl::auxiliary_independent<index_space::edges,
        1 /* dimension */,
        2 /* primary dimension */>>;

    static constexpr size_t auxiliary_colorings =
      std::tuple_size<auxiliary>::value;
    using definition = topo::unstructured_impl::simple_definition;
    using communicator = topo::unstructured_impl::mpi_communicator;
  }; // struct coloring_policy
#endif

  static coloring color(const std::string & fname) {
    io::definition<double> mpas_def(fname.c_str());
    const size_t colors{flecsi::processes()};
    auto [naive, ge, c2v, v2c, c2c] = flecsi::topo::unstructured_impl::make_dcrs(mpas_def, 1);
    auto raw = flecsi::util::parmetis::color(naive, colors);
    auto coloring = flecsi::topo::unstructured_impl::distribute(naive, colors, raw);
#if 0
    auto closure = topo::unstructured_impl::closure<coloring_policy>(
      sd, coloring[0], MPI_COMM_WORLD);

    // FIXME: dummy information so that tests pass
    closure.connectivity_sizes.push_back({10, 10});
    return closure;
#endif
    return {};
  }

  static void initialize(flecsi::data::topology_slot<mesh> & s,
                         const coloring & c) {
    (void)s;
    (void)c;
  } // initialize

}; // struct mpas_mesh

} // namespace poisson
