"""
    Symplectic geometric numerical ODE integrators for Hamiltonian dynamics.
"""

using LinearAlgebra

"""
    step_leapfrog!(system::NBodySystem, dt::Float64)

Standard 2nd-order Velocity Verlet / Leapfrog symplectic integrator step.
Symplectic order: 2.
"""
function step_leapfrog!(system::NBodySystem, dt::Float64)
    n = length(system.bodies)
    a1 = compute_accelerations(system)

    # Half-step velocity kick
    for i in 1:n
        for d in 1:3
            system.bodies[i].velocity[d] += 0.5 * dt * a1[i, d]
        end
    end

    # Full-step position drift
    for i in 1:n
        for d in 1:3
            system.bodies[i].position[d] += dt * system.bodies[i].velocity[d]
        end
    end

    # New acceleration
    a2 = compute_accelerations(system)

    # Second half-step velocity kick
    for i in 1:n
        for d in 1:3
            system.bodies[i].velocity[d] += 0.5 * dt * a2[i, d]
        end
    end
end

"""
    step_yoshida4!(system::NBodySystem, dt::Float64)

Yoshida 4th-order symplectic integrator (Yoshida 1990 / Forest-Ruth).
Preserves phase-space volume and exhibits zero secular energy drift over orbital epochs.
Symplectic order: 4.
"""
function step_yoshida4!(system::NBodySystem, dt::Float64)
    # Yoshida coefficients
    w1 = 1.0 / (2.0 - 2.0^(1.0 / 3.0))
    w0 = - (2.0^(1.0 / 3.0)) / (2.0 - 2.0^(1.0 / 3.0))

    c1 = c4 = w1 * 0.5
    c2 = c3 = (w0 + w1) * 0.5

    d1 = d3 = w1
    d2 = w0

    n = length(system.bodies)

    # Stage 1: Drift c1, Kick d1
    for i in 1:n
        for dim in 1:3
            system.bodies[i].position[dim] += c1 * dt * system.bodies[i].velocity[dim]
        end
    end
    a1 = compute_accelerations(system)
    for i in 1:n
        for dim in 1:3
            system.bodies[i].velocity[dim] += d1 * dt * a1[i, dim]
        end
    end

    # Stage 2: Drift c2, Kick d2
    for i in 1:n
        for dim in 1:3
            system.bodies[i].position[dim] += c2 * dt * system.bodies[i].velocity[dim]
        end
    end
    a2 = compute_accelerations(system)
    for i in 1:n
        for dim in 1:3
            system.bodies[i].velocity[dim] += d2 * dt * a2[i, dim]
        end
    end

    # Stage 3: Drift c3, Kick d3
    for i in 1:n
        for dim in 1:3
            system.bodies[i].position[dim] += c3 * dt * system.bodies[i].velocity[dim]
        end
    end
    a3 = compute_accelerations(system)
    for i in 1:n
        for dim in 1:3
            system.bodies[i].velocity[dim] += d3 * dt * a3[i, dim]
        end
    end

    # Stage 4: Final Drift c4
    for i in 1:n
        for dim in 1:3
            system.bodies[i].position[dim] += c4 * dt * system.bodies[i].velocity[dim]
        end
    end
end

"""
    simulate(system::NBodySystem, t_span::Tuple{Float64, Float64}, dt::Float64; method::Symbol = :yoshida4)::SimulationTrajectory

Simulates Hamiltonian orbital mechanics over `t_span` with step size `dt`.
Supported methods:
- `:yoshida4`: 4th-order symplectic integrator (recommended for long-term orbital stability)
- `:leapfrog`: 2nd-order velocity verlet symplectic integrator
"""
function simulate(
    system::NBodySystem,
    t_span::Tuple{Float64, Float64},
    dt::Float64;
    method::Symbol = :yoshida4,
    save_interval::Int = 1,
)::SimulationTrajectory
    @assert t_span[2] > t_span[1] "Final time must be greater than start time"
    @assert dt > 0.0 "Time step dt must be positive"
    @assert method in (:yoshida4, :leapfrog) "Method must be :yoshida4 or :leapfrog"

    t = t_span[1]
    t_end = t_span[2]
    n = length(system.bodies)

    times = Float64[]
    positions = Matrix{Float64}[]
    velocities = Matrix{Float64}[]
    ke_hist = Float64[]
    pe_hist = Float64[]
    te_hist = Float64[]
    am_hist = Vector{Float64}[]

    step_count = 0

    function record_state!(cur_t)
        push!(times, cur_t)
        pos_mat = zeros(Float64, n, 3)
        vel_mat = zeros(Float64, n, 3)
        for i in 1:n
            pos_mat[i, :] = system.bodies[i].position
            vel_mat[i, :] = system.bodies[i].velocity
        end
        push!(positions, pos_mat)
        push!(velocities, vel_mat)
        ke = kinetic_energy(system)
        pe = potential_energy(system)
        push!(ke_hist, ke)
        push!(pe_hist, pe)
        push!(te_hist, ke + pe)
        push!(am_hist, total_angular_momentum(system))
    end

    record_state!(t)

    while t < t_end
        cur_dt = min(dt, t_end - t)

        if method == :yoshida4
            step_yoshida4!(system, cur_dt)
        else
            step_leapfrog!(system, cur_dt)
        end

        t += cur_dt
        step_count += 1

        if step_count % save_interval == 0 || t >= t_end
            record_state!(t)
        end
    end

    return SimulationTrajectory(
        times,
        positions,
        velocities,
        ke_hist,
        pe_hist,
        te_hist,
        am_hist,
    )
end
