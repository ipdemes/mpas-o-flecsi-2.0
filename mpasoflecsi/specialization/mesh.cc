#include "mesh.hh"


namespace mpas {

using namespace flecsi;

mesh::coloring mesh::color(const std::string & fname)
{
    io::definition<double> mpas_def(fname.c_str());
    const size_t colors{processes()};
    auto [naive, ge, c2v, v2c, c2c] = topo::unstructured_impl::make_dcrs(mpas_def, 1);
    auto raw = util::parmetis::color(naive, colors);
    auto [primaries, p2m, m2p] = topo::unstructured_impl::migrate(naive, colors, raw,
                                                                          c2v, v2c, c2c);
    auto closure = topo::unstructured_impl::closure<coloring_policy>(
      mpas_def, colors, raw, primaries, c2v, v2c, c2c, m2p, p2m);
    {
      auto & curr = closure[process()];
      curr.colors = 1;
      curr.cnx_allocs.push_back(
        {mpas_def.entities_crs(2, 0).indices.size(),
         mpas_def.entities_crs(2, 1).indices.size(),
         mpas_def.entities_crs(2, 2).indices.size()});
      curr.cnx_allocs.push_back({
          mpas_def.entities_crs(1, 0).indices.size(),
          mpas_def.entities_crs(1, 1).indices.size(),
          mpas_def.entities_crs(2, 1).indices.size()});
      curr.cnx_allocs.push_back({
          mpas_def.entities_crs(1, 0).indices.size(),
         mpas_def.entities_crs(2, 0).indices.size()});

      curr.cnx_colorings.push_back(
        {mpas_def.entities_crs(2, 0),
         mpas_def.entities_crs(2, 1),
         mpas_def.entities_crs(2, 2),
         mpas_def.entities_crs(1, 0),
         mpas_def.entities_crs(1, 1)});
    }
    return closure[process()];
}


void mesh::initialize(data::topology_slot<mesh> & s,
                      const mesh::coloring & c)
{
    auto & c2v = s->connect_.get<mesh::cells>().get<mesh::vertices>();
    execute<init_cnx, mpi>(c2v(s), c.cnx_colorings[0][0]);
    auto & v2c = s->connect_.get<mesh::vertices>().get<mesh::cells>();
    execute<transpose_cnx, mpi>(v2c(s), c2v(s));

    auto & c2e = s->connect_.get<mesh::cells>().get<mesh::edges>();
    execute<init_cnx, mpi>(c2e(s), c.cnx_colorings[0][1]);
    auto & e2c = s->connect_.get<mesh::edges>().get<mesh::cells>();
    execute<transpose_cnx, mpi>(e2c(s), c2e(s));

    auto & c2c = s->connect_.get<mesh::cells>().get<mesh::cells>();
    execute<init_cnx, mpi>(c2c(s), c.cnx_colorings[0][2]);

    auto & e2v = s->connect_.get<mesh::edges>().get<mesh::vertices>();
    execute<init_cnx, mpi>(e2v(s), c.cnx_colorings[0][3]);
    auto & v2e = s->connect_.get<mesh::vertices>().get<mesh::edges>();
    execute<transpose_cnx, mpi>(v2e(s), e2v(s));

    auto & e2e = s->connect_.get<mesh::edges>().get<mesh::edges>();
    execute<init_cnx, mpi>(e2e(s), c.cnx_colorings[0][4]);
}


void mesh::init_cnx(field<util::id, data::ragged>::mutator<rw, na> e2e,
                    const topo::unstructured_impl::crs & cnx)
{
    for ( std::size_t c{0}; c < cnx.offsets.size() - 1; ++c) {
      auto start = cnx.offsets[c];
      auto size = cnx.offsets[c+1] - start;

      e2e[c].resize(size);
      for (std::size_t i{0}; i < size; ++i) {
        e2e[c][i] = cnx.indices[start + i];
      }
    }
}


void mesh::transpose_cnx(field<util::id, data::ragged>::mutator<rw, na> v2c,
                         field<util::id, data::ragged>::accessor<ro, na> c2v)
{
  for (std::size_t c{0}; c < c2v.size(); ++c) {
    for (std::size_t v{0}; v < c2v[c].size(); ++v) {
      v2c[c2v[c][v]].push_back(c);
    }
  }
}


}
