import SatelliteDynamics
using LinearAlgebra
using OrdinaryDiffEq
using Plots
plotly()

include("satellite_models.jl")

struct OrbitSimulationParameters
    satellite_model::SatelliteModel
    # scales
    distance_scale::Real
    time_scale::Real
    control_scale::Real
    angular_rate_scale::Real
    # gravity fidelity
    n_gravity::Int
    m_gravity::Int
    # time
    start_epoch::SatelliteDynamics.Epoch # the epoch at which the dynamics integration begins
    # flags
    control_type::Symbol # one of :thrust :drag_ratio :torque or :dipole
    magnetic_model::Symbol
    add_solar_radiation_pressure::Bool
    add_sun_thirdbody::Bool
    add_moon_thirdbody::Bool
end

function OrbitSimulationParameters(
    satellite_model::SatelliteModel;
    distance_scale::Real=1.0,
    time_scale::Real=1.0,
    control_scale::Real=1.0,
    angular_rate_scale::Real=1.0,
    n_gravity::Int=2,
    m_gravity::Int=0,
    start_epoch::SatelliteDynamics.Epoch=SatelliteDynamics.Epoch("2018-01-01"),
    control_type::Symbol=:dipole, # one of :thrust :drag_ratio :torque or :dipole
    magnetic_model::Symbol=:IGRF13,
    add_solar_radiation_pressure::Bool=false,
    add_sun_thirdbody::Bool=false,
    add_moon_thirdbody::Bool=false
)
    return OrbitSimulationParameters(
        satellite_model,
        distance_scale,
        time_scale,
        control_scale,
        angular_rate_scale,
        n_gravity,
        m_gravity,
        start_epoch,
        control_type,
        magnetic_model,
        add_solar_radiation_pressure,
        add_sun_thirdbody,
        add_moon_thirdbody
    )
end

""" cartesian_acceleration_torque(x::Array{<:Real}, params::OrbitSimulationParameters, epc::Epoch)
Compute the accelerations and torque experienced by a satellite
"""
function cartesian_acceleration(x::Array{<:Real}, u::Array{<:Real}, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)

    r = x[1:3] # ECI frame position
    v = x[4:6] # Inertial frame velocity

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = SatelliteDynamics.bias_precession_nutation(epc)
    E = SatelliteDynamics.earth_rotation(epc)
    W = SatelliteDynamics.polar_motion(epc)
    R = W * E * PN

    # eltype(x) makes this forward diff friendly
    a = zeros(eltype(x), 3) # inertial frame acceleration

    # spherical harmonic gravity
    a += SatelliteDynamics.accel_gravity(x, R, params.n_gravity, params.m_gravity)

    # thrust control
    if params.control_type == :thrust
        a += u
    end

    # atmospheric drag
    if params.control_type == :drag_ratio
        drag_ratio = u[1]
    else
        drag_ratio = 0
    end
    a_drag_inertial = drag_ratio_acceleration(x, drag_ratio, epc, params)
    a += a_drag_inertial

    # Sun and Moon
    # third body sun
    if params.add_sun_thirdbody
        r_sun = SatelliteDynamics.sun_position(epc)
        a += SatelliteDynamics.accel_thirdbody_sun(x, r_sun)
    end

    # third body moon
    if params.add_moon_thirdbody
        r_moon = SatelliteDynamics.moon_position(epc)
        a += SatelliteDynamics.accel_thirdbody_moon(x, r_moon)
    end

    return a
end

""" cartesian_keplarian_dynamics(x::Vector{<:Real})
The unperturbed keplarian cartesian orbit dynamics.
"""
function cartesian_keplarian_dynamics(x::Vector{<:Real}, u::Vector{<:Real}, ts::Real)
    r = x[1:3]
    rmag = sqrt(r'r)
    v = x[4:6]

    ṙ = v
    v̇ = -((SatelliteDynamics.GM_EARTH / rmag^3) * r)

    return [ṙ; v̇]
end


""" drag_acceleration_torque(airspeed, epc, params)
Compute the resultant forces and torques on a spacecraft due to drag.
See Markley and Crassidis pg 108

 * `airspeed` is the body-frame velocity of the spacecraft wrt the atmosphere (drag acts opposite to airspeed)
 * `epc` is the current epoch
 * `params::OrbitSimulationParameters` is a struct of simulation and model parameters

All units are SI kg-m-s, return values are in the body-frame
"""
function drag_acceleration_torque(x::Vector{<:Real}, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)
    model = params.satellite_model

    q = x[7:10]

    airspeed_inertial = airspeed_from_state(x)
    airspeed = Q(q)' * airspeed_inertial

    force = zeros(eltype(airspeed), 3)
    torque = zeros(eltype(airspeed), 3)

    ρ = SatelliteDynamics.density_harris_priester(epc, x[1:3])
    CD = model.drag_coefficient
    CoM = model.center_of_mass

    for f in model.faces
        if airspeed'f.normal > 0.0
            f_force = -0.5 * ρ * CD * f.area * (airspeed'f.normal) * airspeed
            force += f_force
            torque += cross(f.position + CoM, f_force)
        end
    end
    accel = force / model.mass
    return accel, torque
end

""" drag_acceleration_torque(airspeed, epc, params)
Compute the resultant forces and torques on a spacecraft due to drag.
See Markley and Crassidis pg 108

 * `airspeed` is the body-frame velocity of the spacecraft wrt the atmosphere (drag acts opposite to airspeed)
 * `epc` is the current epoch
 * `params::OrbitSimulationParameters` is a struct of simulation and model parameters

All units are SI kg-m-s, return values are in the body-frame
"""
function drag_ratio_acceleration(x::Vector{<:Real}, drag_ratio::Real, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)
    @assert 0 <= drag_ratio <= 1.0
    model = params.satellite_model

    v = airspeed_from_state(x)
    v_mag = norm(v)

    ρ = SatelliteDynamics.density_harris_priester(epc, x[1:3])
    C_D = model.drag_coefficient
    mass = model.mass
    A = model.min_drag_area + drag_ratio * (model.max_drag_area - model.min_drag_area)

    F_drag = -0.5 * ρ * C_D * A * v_mag .* v
    a_drag = F_drag / mass

    return a_drag
end

""" solar_acceleration_torque(x::Vector{<:Real}, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)
Compute the acceleration and torque due to solar radiation pressure on each face.
Assumes the SRP center of pressure is the geometric center of each face.
Does not account for self-shielding of faces, assumes constant distance between earth and sun (P_sun is constant).
See Montenbruk and Gill eq 3.73 and Markley and Crassidis eq 3.167
"""
function solar_acceleration_torque(x::Vector{<:Real}, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)
    model = params.satellite_model

    r = x[1:3] # ECI frame position
    v = x[4:6] # ECI frame velocity
    q = x[7:10] # body to inertial quaternion

    r_sun = SatelliteDynamics.sun_position(epc)
    r_sun_body = Q(q)'r_sun
    s = r_sun_body / norm(r_sun_body)

    a = zeros(eltype(x), 3)
    τ = zeros(eltype(x), 3)

    nu = SatelliteDynamics.eclipse_conical([r; v], r_sun)
    P_sun = SatelliteDynamics.P_SUN

    for f in model.faces
        cθ = s'f.normal
        ϵ = f.reflectivity_coefficient_ϵ
        if cθ > 0.0
            force_srp = -nu * P_sun * f.area * cθ * ((1 - ϵ) * s + 2 * ϵ * cθ * f.normal)
            τ += cross(f.position + model.center_of_mass, force_srp)
            a += force_srp / model.mass
        end
    end

    return a, τ
end

""" airspeed_from_state(x)
Compute the inertial-frame airspeed wrt to the atmosphere given the satellite state `x`
"""
function airspeed_from_state(x)
    r = x[1:3] # ECI frame position
    v = x[4:6] # Inertial frame velocity

    ω_atmosphere = [0, 0, SatelliteDynamics.OMEGA_EARTH]
    airspeed_inertial = v + cross(ω_atmosphere, r)

    return airspeed_inertial
end

""" magnetic_torque()

"""
function magnetic_torque(x::Vector{<:Real}, m::Vector{<:Real}, epc::SatelliteDynamics.Epoch, params::OrbitSimulationParameters)
    r = x[1:3]
    B = magnetic_B_vector(r, epc, params)
    τ_mag = cross(m, B)
    return τ_mag
end
function magnetic_B_vector(r, epc, params)
    if params.magnetic_model == :dipole
        return dipole_magnetic_field_B(r)
    elseif params.magnetic_model == :IGRF13
        return IGRF13(r, epc)
    else
        error("Magnetic model $(params.magnetic_model) unrecognized")
    end
end

function satellite_orbit_dynamics!(zdot::Array{<:Real}, z::Array{<:Real}, p, t::Real)
    dscale = p.distance_scale
    tscale = p.time_scale
    uscale = p.control_scale

    # Unpack and unscale state vector
    r = z[1:3] * dscale # ECI frame position
    v = z[4:6] * (dscale / tscale) # Inertial frame velocity
    u = z[7:9] * uscale

    # state vector for cartesian_acceleration_torque()
    x = [r; v]

    epc = sim_time_to_epoch(params, t)
    a = cartesian_acceleration(x, u, epc, params)

    ṙ = v
    v̇ = a

    zdot .= [ṙ / (dscale / tscale)
        v̇ / (dscale / tscale^2)
        zeros(eltype(u), 3)]
end

function sim_time_to_epoch(model, t)
    return model.start_epoch + t
end

function null_controller(x::Vector{<:Real}, t, model::SatelliteModel)
    return zeros(eltype(x), 3)
end

function simulate_satellite_orbit(x0::Array{<:Real}, params::OrbitSimulationParameters, tspan::Tuple{<:Real,<:Real}; controller=null_controller, controller_dt=1.0)

    t0 = tspan[1]
    next_controller_update_time = t0 + controller_dt

    function controls_condition(z, t, integrator)
        params = integrator.p
        tscale = params.time_scale

        return ((t * tscale) - next_controller_update_time) / tscale
    end

    function controls_affect!(integrator)
        params = integrator.p

        dscale = params.distance_scale
        tscale = params.time_scale
        uscale = params.control_scale

        z = integrator.u
        t = integrator.t * tscale

        # Unpack and unscale state vector
        r = z[1:3] * dscale # ECI frame position
        v = z[4:6] * (dscale / tscale) # Inertial frame velocity
        u = z[7:9] * uscale # Body frame magnetic dipole

        x = [r; v]

        u = controller(x, t, params.satellite_model)
        integrator.u[7:9] .= u / uscale
        next_controller_update_time = t + controller_dt
    end

    controls_callback = ContinuousCallback(controls_condition, controls_affect!)

    u0 = controller(x0, t0, params.satellite_model)

    dscale = params.distance_scale
    tscale = params.time_scale
    uscale = params.control_scale

    r0 = x0[1:3] # ECI frame position
    v0 = x0[4:6] # Inertial frame velocity

    z0 = [
        r0 / dscale
        v0 / (dscale / tscale)
        u0 / uscale]

    prob = ODEProblem(satellite_orbit_dynamics!, z0, tspan ./ tscale, params)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, callback=controls_callback, progress=true)

    # unpack time, states, and controls
    thist = tscale .* sol.t
    xhist = [sol[1:3, :] * dscale
        sol[4:6, :] * (dscale / tscale)
    ]
    uhist = sol[7:9, :] * uscale

    return xhist, uhist, thist
end

function plot_position(thist, xhist)
    plot(thist / (60 * 60), xhist[1:3, :]' / 1e3, xlabel="Time (hours)", ylabel="Position (km)", label=["x" "y" "z"], linewidth=1.5)
end

function plot_velocity(thist, xhist)
    plot(thist / (60 * 60), xhist[4:6, :]' / 1e3, xlabel="Time (hours)", ylabel="Velocity (km/s)", label=["vx" "vy" "vz"], linewidth=1.5)
end

function plot_controls(thist, uhist)
    plot(thist / (60 * 60), uhist[1:end, :]', xlabel="Time (hours)", ylabel="Controls (Am²)", label=["mx" "my" "mz"], linewidth=1.5)
end