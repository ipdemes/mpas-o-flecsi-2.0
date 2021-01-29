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


namespace mpas {

struct mpas_mesh : flecsi::topo::specialization<flecsi::topo::unstructured, unstructured> {

}; // struct mpas_mesh

} // namespace poisson
