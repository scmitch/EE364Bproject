# Model for the 2-link nonplanar manipulator
# AA 203 final project, Emiko Soroka
# Modified from astrobee.jl
# Almost all the code is from astrobee.jl except f_dyn, A_dyn, B_dyn.

include("./polygonal_obstacles.jl")

export Manipulator

mutable struct Manipulator

	# State (r, v) and control (u, Γ) dimensions
    x_dim
    u_dim

    # Dynamics and linearized dynamics
    f
    A
    B

    # Model constants
    l1 # Length of first part
    l2 # Length of second part
    model_radius # "Ball" around end effector

    # Problem settings
    dimLinearConstraintsU
    dimSecondOrderConeConstraintsU
    x_init
    x_final
    tf
    h
    xMin
    xMax
    uMin
    uMax
    true_cost_weight

    # Cylindrical obstacles (modeled by a center (x,y) and a radius r) and polygon obstacles (not used in this example)
    obstacles
    poly_obstacles

    # GuSTO parameters
    Delta0
    omega0
    omegamax
    # threshold for constraints satisfaction : constraints <= epsilon
    epsilon
    epsilon_xf_constraint
    rho0
    rho1
    beta_succ
    beta_fail
    gamma_fail
    convergence_threshold # in %
end


function Manipulator()
	x_dim = 5 # [x1, x2, x3, θ₁, θ₂]
	u_dim = 2 # Control θ_dot for θ₁ and θ₂

	# Model parameters
	l1 = 1; l2 = 3;
	model_radius = 0.1;

	dimLinearConstraintsU = 2 # u max/min
	dimSecondOrderConeConstraintsU = 0
	
	# TODO: Be careful changing
  # x[1:3] : Cartesian space coordinates. By inspecting J(q) we get
  # the max/min values x[1:3] can reach.
  # The joint angle max/min is arbitrary.
  xMax = [  l1+l2;   l1+l2;   l2; 2*π; 2*π]
	xMin = [-(l1+l2);-(l1+l2); -l2;   0;   0]
  #ig_limit = 1e6
  #xMax = [ big_limit;   big_limit;   big_limit; 2*π; 2*π]
  #xMin = [-big_limit;-big_limit; -big_limit;   0;   0]

  # arbitrary
	uMax = [1; 1];   uMin = -uMax

	
  # The state is in the augmented space [x, θ] ∈ R⁵
	# TODO: Modify initial x here, in terms of joint angles.
  init_jointspace = [0.; 3.2499] #[π/2;π/2]
  #init_jointspace = [1.75; 3.2499] #[π/2;π/2]
  # 
  x_init = LinearAlgebra.vcat(taskspace(l1, l2, init_jointspace), init_jointspace)
	
  # TODO: modify final x here, in terms of joint angles.
  final_jointspace = [3.033;  0.8666] #[π/4; π/4]
  #final_jointspace = [3.033;  2.8666] #[π/4; π/4]
  x_final = LinearAlgebra.vcat(taskspace(l1, l2, final_jointspace), final_jointspace)


	J0 = get_jacobian(l1, l2, [0;0])
	J0_inv = pinv(J0)

	true_cost_weight = 0.0005

	Delta0 = 100. # Trust region init
  omega0 = 500.
  omegamax = 1.0e6
  epsilon = 1e-3
  epsilon_xf_constraint = 0.
  # TODO: Increase if too many solutions get rejected after you increase N.
  rho0 = 50.0
  rho1 = 100.0
  # Trust region scaling
  beta_succ = 2.
  beta_fail = 0.5
  gamma_fail = 5.
  convergence_threshold = 2.5 # in %

  tf = 1
  # No cylindrical obstacles
  obstacles = []

  # Let's start with one polygonal obstacle
  # TODO: Modify the obstacle here. This one is from Holly's Matlab code.
	sideLengths = [2 2 2];
	corner = [-1 -3 -1];
    center = corner + sideLengths./2
	
	obstacle = PolygonalObstacle(reshape(center, length(center)), reshape(sideLengths, length(sideLengths)))
    poly_obstacles = [obstacle,]

	return Manipulator(
		x_dim, u_dim,
		[], [], [],
		l1, l2,
		model_radius,
		dimLinearConstraintsU,
		dimSecondOrderConeConstraintsU,
		x_init,
		x_final,
		tf, 1, # h, will fix in initialize_trajectory
		xMin,
		xMax,
		uMin,
		uMax,
		true_cost_weight,

		obstacles,
		poly_obstacles,

		Delta0,
		omega0,
		omegamax,

		epsilon,
		epsilon_xf_constraint,
		rho0,
		rho1,
		beta_succ,
		beta_fail,
		gamma_fail,
		convergence_threshold,
		)
end


# From Holly's code
function get_jacobian(l1, l2, q)
	q1 = q[1]
	q2 = q[2]
	J = [-l1*sin(q1)-l2*sin(q1)*sin(q2) l2*cos(q1)*cos(q2)
        l1*cos(q1)+l2*sin(q2)*cos(q1) l2*sin(q1)*cos(q2)
          0                           -l2*sin(q2)]
	return J
end

function taskspace(l1, l2, q)
# This function takes in a joint angle state q and link lengths l and
# computes the end effector position in cartesian space using
# transformation matrices

	x = Array{Float64}(undef, 3)
	
	x[1] = l1*cos(q[1])+l2*sin(q[2])*cos(q[1]);
	x[2] = l1*sin(q[1])+l2*sin(q[1])*sin(q[2]);
	x[3] = l2*cos(q[2]);
	return x
end

# No longer used
#=
# Convert x1, x2, x3 to q1, q2
function inv_taskspace(l1, l2, x)
  q = Array{Float64}(undef, 2)
  # ok so we have x3 = l2 cos(q2)
  # so q2 = arccos(x3/l2)
  # then x2 = sin(q1)*(l1 + l2*sin(q2))
  # so q1 = arcsin(x2/(l1 + l2*sin(q2)))
  q[2] = acos(x[3]/l2)
  q[1] = asin(x[2]/(l1 + l2*sin(q[2])))
  return q
end
=#

# Method that returns the GuSTO parameters (used for set up)

function get_initial_gusto_parameters(m::Manipulator)
    return m.Delta0, m.omega0, m.omegamax, m.epsilon, m.rho0, m.rho1, m.beta_succ, m.beta_fail, m.gamma_fail, m.convergence_threshold
end



# GuSTO is intialized by zero controls, and a straight-line in the state space

function initialize_trajectory(model::Manipulator, N::Int)
  x_dim,  u_dim   = model.x_dim, model.u_dim
  x_init, x_final = model.x_init, model.x_final

  # remember we said we'd fix h? Done.
  model.h = model.tf/N
  
  Q = hcat(range(x_init[4:5], stop=x_final[4:5], length=N)...)
  X = zeros(x_dim, N)
  for i=1:N
    X[:,i] = vcat(taskspace(model.l1, model.l2, Q[:,i]), Q[:,i])
  end
  U = zeros(u_dim, N-1)

  return X, U
end


# Method that returns the convergence ratio between iterations (in percentage)
# The quantities X, U denote the actual solution over time, whereas Xp, Up denote the solution at the previous step over time

function convergence_metric(model::Manipulator, X, U, Xp, Up)
	x_dim = model.x_dim
    N = length(X[1,:])

    # Normalized maximum relative error between iterations
    max_num, max_den = -Inf, -Inf
    for k in 1:N
        val = norm(X[1:x_dim,k] - Xp[1:x_dim,k])
        max_num = val > max_num ? val : max_num

        val = norm(X[1:x_dim,k])
        max_den = val > max_den ? val : max_den
    end

    # Returning percentage error
    return max_num*100.0/max_den
end



# Method that returns the original cost

function true_cost(model::Manipulator, X, U, Xp, Up)
	x_dim, u_dim = model.x_dim, model.u_dim
    cost = 0.

    for k = 1:length(U[1,:])
        cost += sum(U[i,k]^2 for i = 1:u_dim) # This corresponds to ∫ ( || F(t) ||^2 + || T(t) ||^2 ) dt
    end

    return cost
end



# The following methods return the i-th coordinate at the k-th iteration of the various constraints and their linearized versions (when needed)
# These are returned in the form " g(t,x(t),u(t)) <= 0 "

# Method that gathers all the linear control constraints

function control_linear_constraints(model::Manipulator, X, U, Xp, Up, k, i)
    x_dim, u_dim = model.x_dim, model.u_dim
    uMin, uMax = model.uMin, model.uMax

    # Control bounds on u
    if i >= 1 && i<=2
      # uMin <= u
      return uMin[i] - U[i,k]
    elseif i >= 3 && i<=4
      # u <= uMax
      return U[i,k] - uMax[i]
    else
      println("[manipulator.jl::control_linear_constraints] ERROR - too many constraints.")
    end
end

# We don't have any (I think) but I the Astrobee file I started with did
function control_second_order_cone_constraints(model::Manipulator, X, U, Xp, Up, k, i)
    u_dim = model.u_dim
    uMin, uMax = model.uMin, model.uMax

    #if i == 1
    #  return (uMax[1], U[1:2,k])
    #elseif i==2
     # return (uMax[4], U[4:6,k])
    #else
      println("[manipulator.jl::control_second_order_cone_constraints] ERROR - too many constraints.")
    #end
end

# State bounds and trust-region constraints (these are all convex constraints)

function state_max_convex_constraints(model::Manipulator, X, U, Xp, Up, k, i)
    return ( X[i, k] - model.xMax[i] )
end



function state_min_convex_constraints(model::Manipulator, X, U, Xp, Up, k, i)
    return ( model.xMin[i] - X[i, k] )
end



function trust_region_max_constraints(model::Manipulator, X, U, Xp, Up, k, i, Delta)
    return ( (X[i, k] - Xp[i, k]) - Delta )
end



function trust_region_min_constraints(model::Manipulator, X, U, Xp, Up, k, i, Delta)
    return ( -Delta - (X[i, k] - Xp[i, k]) )
end



# Method that checks whether trust-region constraints are satisifed or not (recall that trust-region constraints are penalized)

function is_in_trust_region(model::Manipulator, X, U, Xp, Up, Delta)
    B_is_inside = true

    for k = 1:length(X[1,:])
        for i = 1:model.x_dim
            if trust_region_max_constraints(model, X, U, Xp, Up, k, i, Delta) > 0.
                B_is_inside = false
            end
            if trust_region_min_constraints(model, X, U, Xp, Up, k, i, Delta) > 0.
                B_is_inside = false
            end
        end
    end

    return B_is_inside
end


# Initial and final conditions on state variables

function state_initial_constraints(model::Manipulator, X, U, Xp, Up)
    return ( X[:,1] - model.x_init )
end



function state_final_constraints(model::Manipulator, X, U, Xp, Up)
    return ( X[:,end] - model.x_final )
end






# Methods that return the cylindrical obstacle-avoidance constraint and its lienarized version
# Here, a merely classical distance function is considered

function obstacle_constraint(model::Manipulator, X, U, Xp, Up, k, obs_i,
                                                 obs_type::String="sphere")
    #  obs_type    : Type of obstacles, can be 'sphere' or 'poly'
    #if obs_type=="sphere"
    #  p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
    #  bot_radius        = model.model_radius
    #  total_radius      = obs_radius + bot_radius
    #  p_k  = X[1:3, k]
    #  
    #  dist = norm(p_k - p_obs, 2)
    #  constraint = -( dist - total_radius )
    if obs_type=="poly"
      obs            = model.poly_obstacles[obs_i]
      p_k            = X[1:3, k]
      dist_prev, pos = signed_distance_with_closest_point_on_surface(p_k, obs)
      constraint = -( dist_prev - model.model_radius )
    else
      print("[astrobee.jl::obstacle_constraint_convexified] Unknown obstacle type.")
    end

    return constraint
end



function obstacle_constraint_convexified(model::Manipulator, X, U, Xp, Up, k, obs_i,
                                                 obs_type::String="poly")
    #  obs_type    : Type of obstacles, can be 'sphere' or 'poly'
    #if obs_type=="sphere"
    #  p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
    #  bot_radius        = model.model_radius
    #  total_radius      = obs_radius + bot_radius
    #  p_k, p_kp         = X[1:3, k], Xp[1:3, k]
    #  
    #  dist_prev = norm(p_kp - p_obs, 2)
    #  n_prev    = (p_kp-p_obs) / dist_prev
    #  constraint = -( dist_prev - total_radius + sum(n_prev[i] * (p_k[i]-p_kp[i]) for i=1:3) )

    # Get x,y,z coordinates from θ₁, θ₂
    Xpk = Xp[1:3,k]
    Xk = X[1:3,k]

    if obs_type=="poly"
      obs            = model.poly_obstacles[obs_i]
#      p_k, p_kp      = Xk, Xpk
      dist_prev, pos = signed_distance_with_closest_point_on_surface(Xpk, obs)

      n_prev = (Xpk-obs.c[1:3]) / norm((Xpk-obs.c[1:3]),2)
      constraint = -( dist_prev - model.model_radius + sum(n_prev[i] * (Xk[i]-Xpk[i]) for i=1:3) )
    else
      print("[manipulator.jl::obstacle_constraint_convexified] Unknown obstacle type.")
    end
    
    return constraint
end




# Dynamics and linearization

# The following methods return the dynamical constraints and their linearized versions
# These are returned as time-discretized versions of the constraints " x' - f(x,u) = 0 " or
# " x' - A(t)*(x - xp) - B(t)*(u - up) = 0 ", respectively



# Method that returns the dynamical constraints and their linearized versions all at once

function compute_dynamics(model::Manipulator, Xp, Up)
    N = length(Xp[1,:])

    f_all, A_all, B_all = [], [], []

    for k in 1:N-1
        x_k = Xp[:,k]
        u_k = Up[:,k]

        f_dyn_k, A_dyn_k, B_dyn_k = f_dyn(x_k, u_k, model), A_dyn(x_k, u_k, model), B_dyn(x_k, u_k, model)

        push!(f_all, f_dyn_k)
        push!(A_all, A_dyn_k)
        push!(B_all, B_dyn_k)
    end

    return f_all, A_all, B_all
end


# These methods return the dynamics and its linearizations with respect to the state (matrix A(t)) and the control (matrix B(t)), respectively
# Here is the part we actually modified.

function f_dyn(x::Vector, u::Vector, model::Manipulator)
  x_dim = model.x_dim
  # Nobody should ever write x_p because "p" is ambiguously "plus" and "prev".
  # Of course, in this code we use BOTH.
  #     Xp, Up: previous trajectories
  #     x_p: next discrete-time state xₜ₊₁
  # x_p = zeros(x_dim) # NO! We will resist!
  x_plus = zeros(x_dim) # this is xₜ₊₁

  # We have the state = [x, q]
  # and xₜ₊₁ = xₜ + h*J(q)*q_dot
  #     qₜ₊₁ = qₜ + h*q_dot
  x_plus[1:3] = x[1:3] + model.h.* get_jacobian(model.l1, model.l2, x[4:5])*u
  x_plus[4:5] = x[4:5] + model.h.*u

  return x_plus
end

# How to linearize?
#
# We have x = f(q), x_dot = ∂f/∂q
# we want f(x,u,t) ≈ f(x₀, u₀, t) + ∂f/∂x(x₀,u₀,t)(x-x₀) + ∂f∂u(x₀, u₀, t)(u-u₀)
# and we have u = q_dot

function A_dyn(x::Vector, u::Vector, model::Manipulator)
  x_dim = model.x_dim
  A = zeros(x_dim, x_dim)
  q = x[4:5]
  
  # ∂f/∂x only has dependence on last 2 states: q₁, q₂
  A[1:3,4:5] = model.h*get_jacobian(model.l1, model.l2, q)
  
  # since qₜ₊₁ = qₜ + h*q_dot, partial w/r/t q is I
  A[4:5,4:5] = I(model.u_dim)

  return A
end


function B_dyn(x::Vector, u::Vector, model::Manipulator)
  x_dim, u_dim = model.x_dim, model.u_dim

  B = zeros(x_dim, u_dim)  
  q = x[4:5]

  # xₜ₊₁ = xₜ + h*J(q)*q_dot, so partial q/r/t q_dot is h*J(q)
  B[1:3,:] = model.h*get_jacobian(model.l1, model.l2, q)

  # qₜ₊₁ = qₜ + h*q_dot, so partial w/r/t q_dot is Ih
  B[4:5,:] = I(model.u_dim)*model.h

  return B
end

