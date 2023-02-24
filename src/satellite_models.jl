using LinearAlgebra

struct SatelliteFace
    area::Real # (m^2)
    normal::Vector{<:Real}
    position::Vector{<:Real} # (m) relative to geometric center
    reflectivity_coefficient_ϵ::Real # between 0 and 1, 0 = totally absorbing, 1 = totally reflecting
end

function SatelliteFace(
    area,
    normal,
    position;
    reflectivity_coefficient_ϵ=0.1
)
    return SatelliteFace(
        area,
        normal,
        position,
        reflectivity_coefficient_ϵ
    )
end

struct SatelliteModel
    # fixed parameters
    mass::Real # kg
    inertia::Matrix{<:Real}
    center_of_mass::Vector{<:Real} # (m)
    faces::Vector{SatelliteFace}
    drag_coefficient::Real # C_D, dimensionless
    min_drag_area::Real
    max_drag_area::Real
    max_dipoles::Vector{<:Real} # Am², the maximum dipole that can be produced along each axis
end

dove_faces = [
    SatelliteFace(0.1^2, [0.0, 0.0, 1.0], [0.0, 0.0, 0.3 / 2]) # Top +Z
    SatelliteFace(0.1^2, [0.0, 0.0, -1.0], [0.0, 0.0, -0.3 / 2]) # Bottom -Z
    SatelliteFace(0.1 * 0.3, [0.0, 1.0, 0.0], [0.0, 0.05, 0.0]) # Side +Y
    SatelliteFace(0.1 * 0.3, [0.0, -1.0, 0.0], [0.0, -0.05, 0.0]) # Side -Y
    SatelliteFace(0.1 * 0.3, [1.0, 0.0, 0.0], [0.05, 0.0, 0.0]) # Side +X
    SatelliteFace(0.1 * 0.3, [-1.0, 0.0, 0.0], [-0.05, 0.0, 0.0]) # Side -X
    SatelliteFace(0.3 * 0.3, [1.0, 0.0, 0.0], [0.05, (0.05 + 0.3), 0.0]) # Solar 1 +X face
    SatelliteFace(0.3 * 0.3, [-1.0, 0.0, 0.0], [0.05, (0.05 + 0.3), 0.0]) # Solar 1 -X face
    SatelliteFace(0.3 * 0.3, [1.0, 0.0, 0.0], [0.05, -(0.05 + 0.3), 0.0]) # Solar 2 +X face
    SatelliteFace(0.3 * 0.3, [-1.0, 0.0, 0.0], [0.05, -(0.05 + 0.3), 0.0]) # Solar 2 -X face
]

dove_model = SatelliteModel(
    6.0,
    diagm([1.0, 2.0, 3.0]),
    [0.0, 0.0, 0.0],
    dove_faces,
    2.2,
    0.1^2, # min area
    (2 * 0.3 * 0.3 + 0.1 * 0.3), # max area
    zeros(3),
)

no_drag_model = SatelliteModel(
    1.0,
    diagm([1.0, 1.0, 1.0]),
    [0.0, 0.0, 0.0],
    [],
    0.0,
    0.0,
    0.0,
    zeros(3),
)
