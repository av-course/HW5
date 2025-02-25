using Infiltrator

"""
Unicycle model
"""
function f(x, u, ω, Δ)
    v = x[3]+0.5*Δ*(u[1]+ω[1])
    θ = x[4]+0.5*Δ*(u[2]+ω[2])
    x + Δ * [v*cos(θ), v*sin(θ), u[1]+ω[1], u[2]+ω[2]]
end

"""
Jacobian of f with respect to x, evaluated at x,u,ω,Δ.
"""
function jac_fx(x, u, ω, Δ)
    # TODO
    zeros(4,4)
end


"""
Jacobian of f with respect to u, evaluated at x,u,ω,Δ.
"""
function jac_fu(x, u, ω, Δ)
    # TODO
    zeros(4,2)
end

"""
Jacobian of f with respect to ω, evaluated at x,u,ω,Δ.
"""
function jac_fω(x, u, ω, Δ)
    # TODO 
    zeros(4,2)
end

"""
Non-standard measurement model. Can we extract state estimate from these measurements?
"""
function h(x)
    [atan(x[2], x[1]), 
     -cos(x[4])*x[3]*(x[2]-3*x[1])]
end

"""
Jacobian of h with respect to x, evaluated at x.
"""
function jac_hx(x)
    # TODO
    zeros(2,4)
end


"""
Extended kalman filter implementation.

Assume that the 'true' physical update in the world is given by 

xₖ = f(xₖ₋₁, uₖ, ωₖ, Δ), where Δ is the time difference between times k and k-1.

Here, uₖ is the 'true' controls applied to the system. These controls can be assumed to be a random variable,
with probability distribution given by 𝒩 (mₖ, proc_cov) where mₖ is some IMU-like measurement, and proc_cov is a constant covariance matrix.

ωₖ is assumed to be some random disturbance which affects the system. This could be something like wind. This variable is also presumed to be random,
with probability distribution given by 𝒩 (0, dist_cov).

The process model distribution is then approximated as:

P(xₖ | xₖ₋₁, uₖ) ≈ 𝒩 ( Axₖ₋₁ + Buₖ + L*0 + c, Σ̂ )

where 
A = ∇ₓf(μₖ₋₁, mₖ, 0, Δ),
B = ∇ᵤf(μₖ₋₁, mₖ, 0, Δ),
L = ∇ω f(μₖ₋₁, mₖ, 0, Δ),
c = f(μₖ₋₁, mₖ, 0, Δ) - Aμₖ₋₁ - Bmₖ - L*0

μ̂ = Aμₖ₋₁ + Bmₖ + L*0 + c
  = f(μₖ₋₁, mₖ, 0, Δ)
Σ̂ = A Σₖ₋₁ A' + B proc_cov B' + L dist_cov L'

Further, assume that the 'true' measurement generation in the world is given by

zₖ = h(xₖ) + wₖ,

where wₖ is some additive gaussian noise with probability density function given by

𝒩 (0, meas_var).

The measurement model is then approximated as 

P(zₖ | xₖ) ≈ 𝒩 ( C xₖ + d , meas_var )

where 
C = ∇ₓ h(μ̂), 
d = h(μ̂) - Cμ̂

The extended Kalman filter update equations can be implemented as the following:

Σₖ = (Σ̂⁻¹ + C' (meas_var)⁻¹ C)⁻¹
μₖ = Σₖ ( Σ̂⁻¹ μ̂ + C' (meas_var)⁻¹ (zₖ - d) )

"""
function filter(; μ=zeros(4), Σ=Diagonal([5,5,3,1.0]), x0=zeros(4), num_steps=25, meas_freq=0.5, meas_jitter=0.025, meas_var=Diagonal([0.25,0.25]), proc_cov = Diagonal([0.2, 0.1]), dist_cov=Diagonal([0.3,0.3]), rng=MersenneTwister(5), output=true)
    gt_states = [x0,] # ground truth states that we will try to estimate
    timesteps = []
    u_constant = randn(rng) * [5.0, 0.2]
    μs = [μ,]
    Σs = Matrix{Float64}[Σ,]
    zs = Vector{Float64}[]

    u_prev = zeros(2)
    x_prev = x0

    for k = 1:num_steps
        uₖ = u_constant
        mₖ = uₖ + sqrt(proc_cov) * randn(rng, 2) # Noisy IMU measurement.
        Δ = meas_freq + meas_jitter * (2*rand(rng) - 1)
        ω_true = sqrt(dist_cov) * randn(rng, 2)
        xₖ = f(x_prev, uₖ, ω_true, Δ)
        x_prev = xₖ
        u_prev = uₖ
        zₖ = h(xₖ) + sqrt(meas_var) * randn(rng, 2)

        # TODO : perform update on Σ, μ
        # μ = ...
        # Σ = ...

        push!(μs, μ)
        push!(Σs, Σ)
        push!(zs, zₖ)
        push!(gt_states, xₖ)
        push!(timesteps, Δ)
        if output
            println("Ttimestep ", k, ":")
            println("   Ground truth (x,y): ", xₖ[1:2])
            println("   Estimated (x,y): ", μ[1:2])
            println("   Ground truth v: ", xₖ[3])
            println("   estimated v: ", μ[3])
            println("   Ground truth θ: ", xₖ[4])
            println("   estimated θ: ", μ[4])
            println("   measurement received: ", zₖ)
            println("   Uncertainty measure (det(cov)): ", det(Σ))
        end
    end

    (; μs, Σs)
end
