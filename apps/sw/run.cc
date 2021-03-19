#include "mpasoflecsi/tasks/io.hh"
#include "state.hh"
#include "tasks/solve.hh"
#include "run.hh"

namespace mpas { namespace sw { namespace action {

using namespace flecsi;

static auto output_task = mpas::task::field_output_wrapper<inputs::output_fields>::output;


static void advance_timestep(double dt)
{
  execute<task::init_provis>(m, u[prev](m), h[prev](m), tracers[prev](m),
                             h_edge[prev](m), ke[prev](m), pv_edge[prev](m),
                             u[provis](m), h[provis](m), tracers[provis](m),
                             h_edge[provis](m), ke[provis](m), pv_edge[provis](m));

  execute<task::setup_timestep>(m, u[prev](m), h[prev](m), tracers[prev](m),
                                u[next](m), h[next](m), tracers[next](m));

  std::array<double, 4> rk_weights{{dt / 6, dt / 3, dt / 3, dt / 6}};
  std::array<double, 4> rk_substep_weights{{dt / 2, dt / 2, dt, 0}};

  for(std::size_t rkstep = 0; rkstep < rk_weights.size(); rkstep++) {
    // compute tendencies
    execute<task::compute_tend>(m, areaCell(m), dvEdge(m), dcEdge(m),
                                weightsOnEdge(m), h_edge[provis](m),
                                h[provis](m), bottomDepth(m),
                                u[provis](m), ke[provis](m), pv_edge[provis](m),
                                tend_h(m), tend_u(m));
    execute<task::compute_scalar_tend>(m, areaCell(m), dvEdge(m),
                                       h_edge[provis](m), u[provis](m),
                                       tracers[provis](m), tend_tracers(m));

    // compute next substep state
    if(rkstep < 3) {
      execute<task::compute_substep>(m,
                                     u[provis](m), u[prev](m), tend_u(m),
                                     h[provis](m), h[prev](m), tend_h(m),
                                     tracers[provis](m), tracers[prev](m), tend_tracers(m),
                                     rk_substep_weights[rkstep]);
      execute<task::compute_solve_diagnostics>(m, dvEdge(m), dcEdge(m), areaCell(m),
                                               areaTriangle(m), kiteAreasOnVertex(m),
                                               fVertex(m), h[provis](m), u[provis](m),
                                               h_edge[provis](m), h_vertex[provis](m),
                                               circulation(m), vorticity(m), ke[provis](m),
                                               pv_edge[provis](m), pv_vertex(m));
    }

    // accumulate update (for RK4)
    execute<task::accumulate_update>(m, u[next](m), tend_u(m),
                                     h[next](m), tend_h(m), tracers[next](m), tend_tracers(m),
                                     rk_weights[rkstep]);
  }

  // timestep cleanup: decouple new scalar fields
  execute<task::finalize_timestep>(m, u[prev](m), u[next](m), h[next](m),
                                   tracers[next](m));
  execute<task::compute_solve_diagnostics>(m, dvEdge(m), dcEdge(m),
                                           areaCell(m), areaTriangle(m),
                                           kiteAreasOnVertex(m), fVertex(m),
                                           h[next](m), u[next](m), h_edge[next](m),
                                           h_vertex[next](m), circulation(m), vorticity(m),
                                           ke[next](m), pv_edge[next](m), pv_vertex(m));
  execute<task::timeshift>(m, u[prev](m), h[prev](m), tracers[prev](m), h_edge[prev](m),
                           ke[prev](m), pv_edge[prev](m), u[next](m), h[next](m),
                           tracers[next](m), h_edge[next](m), ke[next](m), pv_edge[next](m));
}


int run()
{
  hid_t mfile;

  auto & state = control::state();

  if (state.output()) {
    execute<mpas::task::mesh_output_init, mpi>(m, areaCell(m), dvEdge(m), dcEdge(m),
                                               xCell(m), yCell(m), zCell(m),
                                               latCell(m), lonCell(m), latVertex(m), lonVertex(m),
                                               meshDensity(m), "test.ncdf", &mfile);
    execute<output_task, mpi>(m, mfile, 0, h[curr](m), bottomDepth(m),
                              vorticity(m), pv_vertex(m), tracers[curr](m));
  }

  execute<task::compute_solve_diagnostics>(m, dvEdge(m), dcEdge(m),
                                           areaCell(m), areaTriangle(m), kiteAreasOnVertex(m),
                                           fVertex(m), h[curr](m), u[curr](m), h_edge[curr](m),
                                           h_vertex[curr](m), circulation(m), vorticity(m),
                                           ke[curr](m), pv_edge[curr](m), pv_vertex(m));

  state.elapsed() = 0;
  double dt = 40;
  std::size_t time_cnt{1};
  for (std::size_t i{1}; i <= state.steps(); ++i) {
    flog(info) << "Running timestep " << i << std::endl;
    advance_timestep(dt);
    if (state.output_step(i)) {
      execute<output_task, mpi>(m, mfile, time_cnt, h[curr](m), bottomDepth(m),
                                vorticity(m), pv_vertex(m), tracers[curr](m));
      ++time_cnt;
    }
    state.elapsed() += dt;
  }

  if (state.output()) {
    execute<mpas::task::mesh_output_finalize, mpi>(mfile);
  }

  return 0;
}

}}}
