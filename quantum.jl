using Plots, LaTeXStrings

# NOTE: Schrodinger's Equation will be abbreviated as SE in most comments.

x = range(-10, 10, length = 5000)

delta_x = x[begin + 1] - x[begin]

function norm(phi) # . means that the function proceeds element-wise
    # abs2 return the square of the absolute value
    norm = sum(abs2, phi) * delta_x
    phi / sqrt(norm)
end

function complex_plot(x, y, probability = true; kwargs... ) # fix this later
    real_comp = real.(y)
    imag_comp = imag.(y)
    # Look at how to handle *args and **kwargs equivalents in julia
    
    plot(xlims = (-2,2))
    plot!(ylims = (-2,2))

    a = plot!(x, real_comp, label = "Re"; kwargs...)
    b = plot!(x, imag_comp, label = "Im"; kwargs...)

    if probability
        p = plot!(x , abs.(y), label = L"\sqrt{P}") 
        # return a, b, p
        return p # Mutates each time, so p holds a and b as well
    else
        # return a, b
        return b
    end
end

# Creates a gaussian wave packet
# This is one of many different wavefunctions that can be used
# pos is position, mom is momentum, sigma is dispersion relation (roughly)
function wave_packet(;pos = 0, mom = 0, sigma = 0.2)
    # Check if any changes for second exp need to be made for complex numbers
    # Look into @. macro
    norm( exp.(-im * mom .* x) .* exp.(-1 .* (x .- pos).^2 ./ sigma^2) )
end

# Great. Now we need to create the framework to time-evole the function

# Assumes boundary of wavefunction in sampling space is 0
# THIS WORKS
function d_dx2(phi, x=x)
    dphi_dx2 = -2 .* phi
    dphi_dx2[begin:end-1] += phi[begin+1:end]
    dphi_dx2[begin+1:end] += phi[begin:end-1]

    dphi_dx2 ./ delta_x
end


# time derivative of SE, where h is the planck constant and m is the mass
# THIS WORKS
function d_dt(phi, h = 1, m = 100, V = 0)
    im * h / (2 * m) .* d_dx2(phi) - im * V .* phi / h
end


# Since we are not assuming any degree of analytical solvability of wavefunction
# We will use a numeric method
# The simplest method of solving differential equations is Euler's method

# FIX BROADCASTING, @. WILLY-NILLY DOESN'T WORK
function euler(phi, dt)
    phi .+ dt .* d_dt(phi)
end

# Upon consulting a numerical analysis book, I found that Euler's method is slow
# compared to more advanced methods, since it assumes the derivative
# at a point is the derivative over a region, requiring smaller steps for accuracy

# As such, I instead will implement a Runge-Kutta integrator of the 4th order

function rk4(phi, dt; kwargs...)
    k1 = d_dt(phi)
    k2 = d_dt(phi .+ dt/2 .* k1)
    k3 = d_dt(phi .+ dt/2 .* k2)
    k4 = d_dt(phi .+ dt .* k3)

    phi .+ dt/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end



# Look at rewriting this in general
# Check how mutators work
function simulate(phi_sim;
                  method = "rk4",
                  V = 0,
                  steps = 100000,
                  dt = 1e-1,
                  condition = nothing,
                  normalize = true,
                  save_every = 100)

    # SEE about taking if method statements out, and replacing with a single
    # check and then assigning function to some integrator_method
    simulation_steps = [copy(phi_sim)]
    
    for i in 1:steps
        if method == "euler"
            phi_sim = euler(phi_sim, dt)
        elseif method == "rk4"
            phi_sim = rk4(phi_sim, dt)
        else
            println("Method ", method, " unknown")
        end

        if !isnothing(condition)
            phi_sim = condition(phi_sim)
        end
        if normalize
            phi_sim = norm(phi_sim)
        end
        if !isnothing(save_every) && i % save_every == 0
            push!(simulation_steps, copy(phi_sim))
        end
    end

    return simulation_steps
end

# Now to create a scenario. In this case, it *is* easiest to create
# a free particle, as opposed to a infinite well
sim_free = simulate(wave_packet(), dt = 1e-1)


# Now we need to develop consierations for nonzero potentials.


# Now to create a function that animates these situations

function animate(simulation_steps, init_func = nothing)
    anim = @animate for i in 1:100000
        simulation_steps[i]
end

animate(sim_free)

@gif for i in 1:1000
    complex_plot(x, sim_free[i])
end

anim = @animate for i in 1:100
    complex_plot(x, sim_free[10i])
end