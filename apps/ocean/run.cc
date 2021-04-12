#include "state.hh"
#include "tasks/solve.hh"
#include "run.hh"

namespace mpas { namespace ocean { namespace action {

using namespace flecsi;


static void rk4_timestep(double dt)
{
  execute<task::copy_tracers>(m, tracers[curr](m), tracers[next](m), tracers[prev](m));

  std::array<double, 4> rk_weights{{dt / 6, dt / 3, dt / 3, dt / 6}};
  std::array<double, 4> rk_substep_weights{{dt / 2, dt / 2, dt, 0}};

  for(std::size_t rkstep = 0; rkstep < rk_weights.size(); rkstep++) {
//    flog(info) << "Substep " << rkstep + 1 << std::endl;
    execute<task::reinitialize_tracer_tendencies>(m, tracers_tend(m));
    if(control::state().monoAdvect_config()) {
      execute<task::calc_normal_thickness_flux>(m, layerThicknessEdge(m), normalVelocity[curr](m),
                                                normalThicknessFlux(m));
      execute<task::calc_provisional_layer_thicknesses>(m, maxLevelCell(m), areaCell(m), dvEdge(m), edgeSignOnCell(m),
                                                        layerThickness[curr](m), normalThicknessFlux(m), w(m),
                                                        layerThickness[prov](m), hProvInv(m), hNewInv(m), dt);
      execute<task::calc_horiz_tracer_advection_flux_mono_low_order>(m, dvEdge(m), normalThicknessFlux(m),
                                                                     tracers[prev](m), lowOrderFlux(m));
      execute<task::calc_horiz_tracer_advection_flux_mono_high_order>(m, highOrderAdvMask(m), dvEdge(m),
              normalThicknessFlux(m), advCoef(m), advCoef_3rd(m), maxLevelCell(m), nAdvCellsForEdge(m),
              advCellsForEdge(m), tracers[prev](m), highOrderFlux(m));
      execute<task::correct_and_rescale_high_order_mono_horiz_flux>(m, edgeSignOnCell(m), hProvInv(m),
              layerThickness[curr](m), normalThicknessFlux(m), areaCell(m), dvEdge(m), highOrderAdvMask(m),
              maxLevelCell(m), maxLevelEdgeTop(m), lowOrderFlux(m), tracers[prev](m), highOrderFlux(m),
              flxIn(m), flxOut(m), tracers_tend(m), tracerMax(m), tracerMin(m), dt);
      execute<task::accumulate_horiz_tracer_advection_tend_mono_FCT>(m, edgeSignOnCell(m), hProvInv(m),
              layerThickness[curr](m), areaCell(m), maxLevelCell(m), maxLevelEdgeTop(m), highOrderFlux(m),
              tracers[prev](m), tracers[prov](m), tracers_tend(m), dt);
    }
    if(control::state().hmix_config()) {
      execute<task::tracer_hmix_del2_tend>(m, areaCell(m), dcEdge(m), dvEdge(m), meshScalingDel2(m), maxLevelEdgeTop(m),
              edgeSignOnCell(m), layerThicknessEdge(m), tracers[prev](m), tracers_tend(m), control::state().hmix_coef());
    }
    execute<task::step_tracers_forward>(m, tracers_tend(m), tracers[curr](m), tracers[next](m), tracers[prev](m),
                                        rk_weights[rkstep], rk_substep_weights[rkstep]);

  }

  execute<task::update_tracers>(m, tracers[next](m), tracers[curr](m));

}

static void advance_timestep(double dt)
{
  if(control::state().rk4_config()) {
    rk4_timestep(dt);
  } else if(control::state().se_config()) {

  }
}

int run()
{
  auto & state = control::state();

  state.elapsed() = 0;

  double dt = inputs::delta_t.value();
  //double dt = 96.;  

  for (std::size_t i = 1; i <= inputs::nsteps.value(); ++i) {
    flog(info) << "Running timestep " << i << std::endl;
    advance_timestep(dt);
    state.elapsed() += dt;
  }

  return 0;
}

}}}
