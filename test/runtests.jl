using Test
using LinearAlgebra
using HyperionOrbit

@testset "HyperionOrbit.jl Test Suite" begin

    @testset "Center of Mass and Barycentric Frame" begin
        b1 = Body3D("Sun", 2.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        b2 = Body3D("Planet", 0.5, [4.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        sys = NBodySystem([b1, b2])

        r_cm, v_cm = center_of_mass(sys)
        @test r_cm ≈ [0.8, 0.0, 0.0] atol=1e-12
        @test v_cm ≈ [0.0, 0.2, 0.0] atol=1e-12

        drift_to_com_frame!(sys)
        r_cm_after, v_cm_after = center_of_mass(sys)
        @test norm(r_cm_after) < 1e-12
        @test norm(v_cm_after) < 1e-12

        P = total_linear_momentum(sys)
        @test norm(P) < 1e-12
    end

    @testset "Hamiltonian Energy & Symplectic Conservation (Yoshida 4)" begin
        # 2-Body circular orbit: a = 1.0, M1 = 1.0, M2 = 0.01, G = 1.0
        sys = preset_kepler_two_body(mass1=1.0, mass2=0.01, semi_major_axis=1.0, eccentricity=0.1, G=1.0)
        e0 = total_energy(sys)
        l0 = total_angular_momentum(sys)

        # Orbital period T = 2*pi * sqrt(a^3 / G(M1+M2))
        period = 2.0 * pi * sqrt(1.0^3 / (1.0 * 1.01))

        # Simulate 10 full orbits with Yoshida 4th order
        dt = 0.01
        traj = simulate(sys, (0.0, 10.0 * period), dt; method=:yoshida4)

        @test length(traj.times) > 100

        # Maximum relative energy fluctuation
        max_de = maximum(abs.(traj.total_energies .- e0) ./ abs(e0))
        @test max_de < 1e-5

        # Angular momentum conservation (exact to machine precision)
        for L in traj.angular_momenta
            @test norm(L .- l0) / norm(l0) < 1e-10
        end
    end

    @testset "Kepler's Third Law Periodicity" begin
        a = 2.0
        G = 1.0
        m1 = 1.0
        m2 = 1e-5
        sys = preset_kepler_two_body(mass1=m1, mass2=m2, semi_major_axis=a, eccentricity=0.0, G=G)

        # Expected period: T = 2*pi * sqrt(a^3 / G(M1+M2))
        T_expected = 2.0 * pi * sqrt(a^3 / (G * (m1 + m2)))

        init_pos2 = copy(sys.bodies[2].position)

        # Integrate for exactly one orbital period
        dt = 0.005
        traj = simulate(sys, (0.0, T_expected), dt; method=:yoshida4)

        final_pos2 = traj.positions[end][2, :]
        diff = norm(final_pos2 .- init_pos2)
        @test diff < 5e-3
    end

    @testset "Figure-8 Choreography" begin
        sys = preset_figure8_three_body()
        e0 = total_energy(sys)
        l0 = total_angular_momentum(sys)

        # Figure-8 has zero angular momentum due to symmetric choreography
        @test norm(l0) < 1e-6

        # Simulate for 1 full period (T ≈ 6.3259)
        traj = simulate(sys, (0.0, 6.3259), 0.002; method=:yoshida4)

        max_de = maximum(abs.(traj.total_energies .- e0) ./ abs(e0))
        @test max_de < 1e-6
    end

    @testset "Lagrange Points Equilateral Triangle Geometry" begin
        mu = 0.05
        pts = lagrange_points_cr3bp(mu)

        @test haskey(pts, :L1)
        @test haskey(pts, :L2)
        @test haskey(pts, :L3)
        @test haskey(pts, :L4)
        @test haskey(pts, :L5)

        # Primary at (-mu, 0), Secondary at (1 - mu, 0)
        # In normalized CR3BP, distance between primaries is 1.0
        r_p1 = [-mu, 0.0, 0.0]
        r_p2 = [1.0 - mu, 0.0, 0.0]

        d_p1_p2 = norm(r_p2 .- r_p1)
        @test d_p1_p2 ≈ 1.0

        # L4 forms an equilateral triangle: distance to both primaries is exactly 1.0
        d_p1_L4 = norm(pts[:L4] .- r_p1)
        d_p2_L4 = norm(pts[:L4] .- r_p2)

        @test d_p1_L4 ≈ 1.0 atol=1e-12
        @test d_p2_L4 ≈ 1.0 atol=1e-12

        # Collinear points are on x-axis (y = 0, z = 0)
        @test pts[:L1][2] == 0.0
        @test pts[:L2][2] == 0.0
        @test pts[:L3][2] == 0.0
        @test pts[:L1][1] < pts[:L2][1]
    end
end
