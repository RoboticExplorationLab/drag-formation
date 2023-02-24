using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
import SatelliteDynamics

include("../src/satellite_simulator.jl")
include("../src/satellite_models.jl")

x_osc_0 = [400e3 + SatelliteDynamics.R_EARTH, 0.0, deg2rad(98), deg2rad(-1.0), 0.0, 0.0] # a, e, i, Ω, ω, M

params = OrbitSimulationParameters(dove_model;
    distance_scale=1e7,
    time_scale=1e4,
    control_scale=10,
    n_gravity=2,
    m_gravity=0,
    add_sun_thirdbody=true,
    add_moon_thirdbody=true)

tspan = (0.0, 60 * 60.0)
x0 = SatelliteDynamics.sOSCtoCART(x_osc_0; use_degrees=false)

xhist, uhist, thist = simulate_satellite_orbit(x0[1:6], params, tspan; controller=null_controller, controller_dt=60.0)

plot_position(thist, xhist)
plot_velocity(thist, xhist)
plot_controls(thist, uhist)