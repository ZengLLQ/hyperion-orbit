"""
    Core data structures for N-Body Hamiltonian orbital dynamics.
"""

using LinearAlgebra

"""
    Body3D

Represents a single celestial body with mass and phase space coordinates in 3D.
"""
mutable struct Body3D
    name::String
    mass::Float64
    position::Vector{Float64} # [x, y, z] in meters or AU
    velocity::Vector{Float64} # [vx, vy, vz] in m/s or AU/day

    function Body3D(name::String, mass::Float64, position::Vector{Float64}, velocity::Vector{Float64})
        @assert length(position) == 3 "Position must be 3-dimensional"
        @assert length(velocity) == 3 "Velocity must be 3-dimensional"
        @assert mass >= 0.0 "Mass must be non-negative"
        new(name, mass, copy(position), copy(velocity))
    end
end

"""
    NBodySystem

Collection of gravitational bodies governed by Newton's law of universal gravitation
with optional Plummer softening to avoid close-encounter singularity divergence.
"""
struct NBodySystem
    bodies::Vector{Body3D}
    G::Float64
    softening_epsilon::Float64

    function NBodySystem(bodies::Vector{Body3D}; G::Float64 = 1.0, softening_epsilon::Float64 = 0.0)
        @assert G > 0.0 "Gravitational constant G must be positive"
        @assert softening_epsilon >= 0.0 "Softening parameter epsilon must be non-negative"
        new(bodies, G, softening_epsilon)
    end
end

"""
    SimulationTrajectory

Recorded time-series trajectory and conserved dynamical invariants over a simulation run.
"""
struct SimulationTrajectory
    times::Vector{Float64}
    # positions[step][body_idx, :]
    positions::Vector{Matrix{Float64}}
    velocities::Vector{Matrix{Float64}}
    kinetic_energies::Vector{Float64}
    potential_energies::Vector{Float64}
    total_energies::Vector{Float64}
    angular_momenta::Vector{Vector{Float64}}
end
