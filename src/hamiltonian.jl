"""
    Hamiltonian mechanics, conservation laws, and gravitational acceleration kernels.
"""

using LinearAlgebra

"""
    compute_accelerations(system::NBodySystem)::Matrix{Float64}

Computes the gravitational acceleration vector for every body in the system
using Plummer softening:
    a_i = sum_{j != i} G * m_j * (r_j - r_i) / (||r_j - r_i||^2 + eps^2)^(3/2)
Returns an (N x 3) matrix of accelerations.
"""
function compute_accelerations(system::NBodySystem)::Matrix{Float64}
    n = length(system.bodies)
    acc = zeros(Float64, n, 3)
    eps2 = system.softening_epsilon^2

    for i in 1:n
        ri = system.bodies[i].position
        for j in 1:n
            if i != j
                rj = system.bodies[j].position
                mj = system.bodies[j].mass
                dx = rj[1] - ri[1]
                dy = rj[2] - ri[2]
                dz = rj[3] - ri[3]
                dist2 = dx^2 + dy^2 + dz^2 + eps2
                dist = sqrt(dist2)
                inv_dist3 = 1.0 / (dist * dist2)
                factor = system.G * mj * inv_dist3

                acc[i, 1] += factor * dx
                acc[i, 2] += factor * dy
                acc[i, 3] += factor * dz
            end
        end
    end
    return acc
end

"""
    kinetic_energy(system::NBodySystem)::Float64

Total kinetic energy: T = 0.5 * sum(m_i * ||v_i||^2).
"""
function kinetic_energy(system::NBodySystem)::Float64
    ke = 0.0
    for b in system.bodies
        v2 = b.velocity[1]^2 + b.velocity[2]^2 + b.velocity[3]^2
        ke += 0.5 * b.mass * v2
    end
    return ke
end

"""
    potential_energy(system::NBodySystem)::Float64

Total gravitational potential energy with Plummer softening:
U = - sum_{i < j} G * m_i * m_j / sqrt(||r_i - r_j||^2 + eps^2).
"""
function potential_energy(system::NBodySystem)::Float64
    pe = 0.0
    n = length(system.bodies)
    eps2 = system.softening_epsilon^2

    for i in 1:(n - 1)
        ri = system.bodies[i].position
        mi = system.bodies[i].mass
        for j in (i + 1):n
            rj = system.bodies[j].position
            mj = system.bodies[j].mass
            dx = rj[1] - ri[1]
            dy = rj[2] - ri[2]
            dz = rj[3] - ri[3]
            dist = sqrt(dx^2 + dy^2 + dz^2 + eps2)
            pe -= (system.G * mi * mj) / dist
        end
    end
    return pe
end

"""
    total_energy(system::NBodySystem)::Float64

Total Hamiltonian energy H = T + U.
"""
function total_energy(system::NBodySystem)::Float64
    return kinetic_energy(system) + potential_energy(system)
end

"""
    total_angular_momentum(system::NBodySystem)::Vector{Float64}

Total angular momentum vector: L = sum(m_i * (r_i x v_i)).
"""
function total_angular_momentum(system::NBodySystem)::Vector{Float64}
    L = zeros(Float64, 3)
    for b in system.bodies
        # cross product r x v
        lx = b.position[2] * b.velocity[3] - b.position[3] * b.velocity[2]
        ly = b.position[3] * b.velocity[1] - b.position[1] * b.velocity[3]
        lz = b.position[1] * b.velocity[2] - b.position[2] * b.velocity[1]
        L[1] += b.mass * lx
        L[2] += b.mass * ly
        L[3] += b.mass * lz
    end
    return L
end

"""
    total_linear_momentum(system::NBodySystem)::Vector{Float64}

Total linear momentum vector: P = sum(m_i * v_i).
"""
function total_linear_momentum(system::NBodySystem)::Vector{Float64}
    P = zeros(Float64, 3)
    for b in system.bodies
        P[1] += b.mass * b.velocity[1]
        P[2] += b.mass * b.velocity[2]
        P[3] += b.mass * b.velocity[3]
    end
    return P
end

"""
    center_of_mass(system::NBodySystem)::Tuple{Vector{Float64}, Vector{Float64}}

Calculates barycenter position and velocity vectors.
"""
function center_of_mass(system::NBodySystem)::Tuple{Vector{Float64}, Vector{Float64}}
    m_total = sum(b.mass for b in system.bodies)
    @assert m_total > 0.0 "Total mass must be positive"

    r_cm = zeros(Float64, 3)
    v_cm = zeros(Float64, 3)

    for b in system.bodies
        r_cm .+= b.mass .* b.position
        v_cm .+= b.mass .* b.velocity
    end

    return (r_cm ./ m_total, v_cm ./ m_total)
end

"""
    drift_to_com_frame!(system::NBodySystem)

Shifts coordinate frame to the system center of mass and removes net linear drift.
"""
function drift_to_com_frame!(system::NBodySystem)
    r_cm, v_cm = center_of_mass(system)
    for b in system.bodies
        b.position .-= r_cm
        b.velocity .-= v_cm
    end
end
