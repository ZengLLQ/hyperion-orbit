"""
    Astronomical and mathematical celestial mechanics presets.
"""

using LinearAlgebra

"""
    preset_kepler_two_body(; mass1 = 1.0, mass2 = 0.001, semi_major_axis = 1.0, eccentricity = 0.1, G = 1.0)::NBodySystem

Constructs a two-body Keplerian elliptic orbit system.
At periapsis:
    r_peri = a * (1 - e)
    v_peri = sqrt(G * (m1 + m2) * (1 + e) / (a * (1 - e)))
"""
function preset_kepler_two_body(;
    mass1::Float64 = 1.0,
    mass2::Float64 = 0.001,
    semi_major_axis::Float64 = 1.0,
    eccentricity::Float64 = 0.0,
    G::Float64 = 1.0,
)::NBodySystem
    r_peri = semi_major_axis * (1.0 - eccentricity)
    v_peri = sqrt(G * (mass1 + mass2) * (1.0 + eccentricity) / r_peri)

    b1 = Body3D("Primary", mass1, [0.0, 0.0, 0.0], [0.0, -mass2 * v_peri / mass1, 0.0])
    b2 = Body3D("Secondary", mass2, [r_peri, 0.0, 0.0], [0.0, v_peri, 0.0])

    sys = NBodySystem([b1, b2], G = G)
    drift_to_com_frame!(sys)
    return sys
end

"""
    preset_figure8_three_body()::NBodySystem

Periodic choreography of 3 equal masses moving on a figure-8 trajectory
discovered by Chenciner & Montgomery (2000).
G = 1.0, m1 = m2 = m3 = 1.0.
"""
function preset_figure8_three_body()::NBodySystem
    # Initial conditions (Simó 2002 precision)
    x1 = 0.97000436
    y1 = -0.24308753
    vx3 = -0.93240737
    vy3 = -0.86473146

    b1 = Body3D("Body1", 1.0, [x1, y1, 0.0], [-0.5 * vx3, -0.5 * vy3, 0.0])
    b2 = Body3D("Body2", 1.0, [-x1, -y1, 0.0], [-0.5 * vx3, -0.5 * vy3, 0.0])
    b3 = Body3D("Body3", 1.0, [0.0, 0.0, 0.0], [vx3, vy3, 0.0])

    return NBodySystem([b1, b2, b3], G = 1.0, softening_epsilon = 0.0)
end

"""
    preset_sun_earth_jupiter()::NBodySystem

Solar system 3-body configuration with Sun, Earth, and Jupiter in normalized astronomical units.
G = 4*pi^2 (AU^3 / (yr^2 * M_sun)).
"""
function preset_sun_earth_jupiter()::NBodySystem
    G_astro = 4.0 * pi^2

    m_sun = 1.0
    m_earth = 3.003e-6
    m_jupiter = 9.543e-4

    r_earth = 1.0 # 1 AU
    v_earth = sqrt(G_astro * m_sun / r_earth)

    r_jupiter = 5.204 # 5.2 AU
    v_jupiter = sqrt(G_astro * m_sun / r_jupiter)

    sun = Body3D("Sun", m_sun, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    earth = Body3D("Earth", m_earth, [r_earth, 0.0, 0.0], [0.0, v_earth, 0.0])
    jupiter = Body3D("Jupiter", m_jupiter, [0.0, r_jupiter, 0.0], [-v_jupiter, 0.0, 0.0])

    sys = NBodySystem([sun, earth, jupiter], G = G_astro)
    drift_to_com_frame!(sys)
    return sys
end

"""
    lagrange_points_cr3bp(mu::Float64)::Dict{Symbol, Vector{Float64}}

Calculates coordinates of the 5 Lagrange equilibrium points in normalized
Circular Restricted Three-Body Problem (CR3BP) rotating frame.
`mu = m2 / (m1 + m2)`.
"""
function lagrange_points_cr3bp(mu::Float64)::Dict{Symbol, Vector{Float64}}
    @assert 0.0 < mu < 0.5 "Mass ratio mu must be in (0, 0.5)"

    # L4 and L5 form equilateral triangles with the primaries
    x_L45 = 0.5 - mu
    y_L4 = sqrt(3.0) / 2.0
    y_L5 = -sqrt(3.0) / 2.0

    # L1, L2, L3 are roots of Euler collinear polynomial
    # For small mu (Hill sphere approx):
    gamma = (mu / 3.0)^(1.0 / 3.0)
    x_L1 = 1.0 - mu - gamma + (1.0 / 3.0) * gamma^2 + (1.0 / 9.0) * gamma^3
    x_L2 = 1.0 - mu + gamma + (1.0 / 3.0) * gamma^2 - (1.0 / 9.0) * gamma^3
    x_L3 = -(1.0 + (5.0 / 12.0) * mu)

    return Dict(
        :L1 => [x_L1, 0.0, 0.0],
        :L2 => [x_L2, 0.0, 0.0],
        :L3 => [x_L3, 0.0, 0.0],
        :L4 => [x_L45, y_L4, 0.0],
        :L5 => [x_L45, y_L5, 0.0],
    )
end
