# Literally just plt_solutions from astrobee_plotting.jl but
# modified

using Plots

# Python plotting with matplotlib
using PyCall, LaTeXStrings
# using PyPlot
import PyPlot; const plt = PyPlot


function plot_solutions(scp_problem::GuSTOProblem, model, X_all, U_all)
    N = length(X_all)

    idx = [1,2]
    local fig
    fig = plot(X_all[1][idx[1],:], X_all[1][idx[2],:],
        label="iter = 0",linewidth=2,
#        xlims=(-0.5,6.),ylims=(-0.5,6.5),
        xlabel="x",ylabel="y",legend=true)
    
    for iter = 2:length(X_all)
        X = X_all[iter]
        plot!(fig, X[idx[1],:], X[idx[2],:],
            label="iter = $(iter - 1)",linewidth=2)
    end

    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plot_circle(p_obs[idx], obs_radius; color=:red, fig=fig)
    end

    return fig
end


function plt_3D_solutions(scp_problem::GuSTOProblem, model, X_all, U_all;
                        xlims=[-0.5,3.], ylims=[0.0,6.], figsize=(8,6), B_plot_labels=true)
    N = length(X_all)

    font = Plots.font("Helvetica", 16)
    titlefont = Plots.font("Helvetica", 28)
    smallfont = Plots.font("Helvetica", 1)

    pyplot(grid=true, guidefont=font, titlefont=titlefont,
        xtickfont=font, ytickfont=font, ztickfont=font, leg=false)
        # xtickfont=smallfont, ytickfont=smallfont, ztickfont=smallfont, leg=false)

    idx = [1,2,3]
    camera = (40, 40)
    idx = [1,2,3]

    

    xticks = collect(-4:1:4)
    yticks = collect(-4:1:4)
    zticks = collect(-4:1:4)

    local fig
    # Main trajectory plot
    # Plot SCP solutions
    fig = plot(X_all[1][idx[1],:], X_all[1][idx[2],:],X_all[1][idx[3],:],
        linewidth=2, label="Initializer", camera=camera,
        xlims=(-4,4), ylims=(-4,4),zlims=(-4,4),
        xticks=xticks,yticks=yticks,zticks=zticks,
        #title="Manipulator SCP Trajectories",

        )
    ax = plt.gca()
    
    for iter = 2:length(X_all)
        X = X_all[iter]
        fig = plot!(fig, X[idx[1],:], X[idx[2],:], X[idx[3],:],
                    label="Iterate $(iter - 1)", linewidth=2,legend=true)
    end

    ## ---------------------- ##
    ## ----- ENVIRONMENT ---- ##

    ## ----- obstacles ------ ##
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plt_circle(ax, p_obs[idx], obs_radius; color="r", alpha=0.4)
    end
    plt.plot(Inf*[1,1],Inf*[1,1], "r-") # for legend

    
    # Plot the poly obstacles
    for obs in model.poly_obstacles
        center, widths = obs.c, 2. * Vector([obs.dx,obs.dy,obs.dz])
        
       fig = plt_rectangle(fig, center[idx], widths[idx], color="red", alpha=0.4)
    end
    # Plot the torus
    fig = plot_torus(fig)


    # Settings / Style / Parameters
    PyPlot.rc("text", usetex=true)
    rcParams = PyDict(plt.matplotlib["rcParams"])
    rcParams["font.size"] = 22
    #rcParams["font.family"] = "Helvetica"
    # plt.xlim(xlims)
    # plt.ylim(ylims)

    #if B_plot_labels
        #plt.title("Manipulator SCP Trajectories")
        #ax.title.set_position([.5, 1.01])
        #plt.xlabel("X", fontsize=26)
        #plt.ylabel("Y", fontsize=26, rotation="horizontal")    
        #plt.legend(loc="left", labelspacing=0.1)
        #plt.legend()
    #end
    #plt.grid(alpha=0.3)

    ax.legend()
    return fig
end



function plot3D_final_solution(scp_problem::GuSTOProblem, model, X, U)
    N = length(X_all)

    font = Plots.font("Helvetica", 16)
    titlefont = Plots.font("Helvetica", 28)
    smallfont = Plots.font("Helvetica", 1)

    pyplot(grid=true, guidefont=font, titlefont=titlefont,
        xtickfont=font, ytickfont=font, ztickfont=font, leg=false)
        # xtickfont=smallfont, ytickfont=smallfont, ztickfont=smallfont, leg=false)

    idx = [1,2,3]
    camera = (40, 40)

    xlims = (-4, 4)
    ylims = (-4, 4)
    zlims = (-4, 4)

    X = copy(X)
    idx_in =            X[idx[1],:].<= xlims[2]
    idx_in = idx_in .& (X[idx[1],:].>= xlims[1])
    idx_in = idx_in .& (X[idx[2],:].<= ylims[2])
    idx_in = idx_in .& (X[idx[2],:].>= ylims[1])
    X = X[:, idx_in]

    xticks = collect(-4:1:4)
    yticks = collect(round(ylims[1],digits=1):1.:round(ylims[2],digits=1))
    zticks = collect(-4:1:round(zlims[2],digits=1))

    local fig
    # Main trajectory plot
    
    fig = plot(X[idx[1],:],X[idx[2],:],X[idx[3],:],
        linewidth=2,
        #title="Final Trajectory",
        color=:blue,
        xlims=xlims, ylims=ylims,zlims=zlims,
        xticks=xticks,yticks=yticks,zticks=zticks,
        # xlabel="X", ylabel="Y", zlabel="Z",
        camera = camera,
        fontfamily="Helvetica",margin=5Plots.mm,left_margin = 8Plots.mm)

    ax = plt.gca()
    
    # Spherical obstacles
    # color = :red
    color = RGB(255. / 255,  120. / 255,   120. / 255)
    for obs_i = 1:length(model.obstacles)
        p_obs, obs_radius = model.obstacles[obs_i][1], model.obstacles[obs_i][2]
        plot_sphere(p_obs[idx], obs_radius; color=color, fig=fig)
    end
    for obs_i = 1:length(model.obstacles)
        # sphere projection
        p_obs_2d, obs_rad = model.obstacles[obs_i][1][1:2], model.obstacles[obs_i][2]
        u = LinRange(-pi/2.0,pi/2.0,50)
        v = LinRange(0,2*pi,100)
        x, y, z = [], [], []
        for i = 1:length(u)
            for j = 1:length(v)
                push!(x,p_obs_2d[1] + obs_rad*cos(u[i])*cos(v[j]))
                push!(y,p_obs_2d[2] + obs_rad*cos(u[i])*sin(v[j]))
                push!(z,zlims[1])
            end
        end
        fig = plot!(fig,x,y,z,linetype=:surface,colorbar=false,
            seriestype = [:shape,],
            c  = color,alpha = 0.2)
    end

    
    for obs in model.poly_obstacles
        center, widths = obs.c, 2. * Vector([obs.dx,obs.dy,obs.dz])
        
       fig = plt_rectangle(fig, center[idx], widths[idx], color="red", alpha=0.4)
    end

    # Plot the torus
    fig = plot_torus(fig)

    # Trajectory
    fig = plot!(X[idx[1],:],X[idx[2],:],zlims[1]*ones(length(X[idx[3],:])),
        linewidth=2,color=:gray,alpha=0.7)
    fig = scatter!(X[idx[1],:],X[idx[2],:],zlims[1]*ones(length(X[idx[3],:])),
        linewidth=4,color=:gray,markerstrokecolor=:darkgray,
        markersize=6, alpha=0.7)
    fig = scatter!(X[idx[1],:],X[idx[2],:],X[idx[3],:],
        linewidth=4,color=:blue,markerstrokecolor=:darkblue,
        markersize=6, alpha=1.)

    xlabel!("x")
    ylabel!("y")
    # for iter = 2:length(X_all)
    #     X = X_all[iter]
    #     plot!(fig,X[idx[1],:],X[idx[2],:],X[idx[3],:],
    #         linewidth=2,label="iter = $(iter - 1)")
    # end

    return fig
end





# No! Don't look here!
# There isn't a plot cube function in Plots.jl
# If you look past here, only pain and sorrow awaits.
# Just scroll back the way you came and pretend this file ends at line 232.

function plt_rectangle(fig, center, widths, additional_w=0; 
                           color="b", alpha=0.1, noFaceColor=false, 
                           label="None")
    """
    Plots a rectangle with a given center and total widths
    arguments:  - center    - (2,) center
                - widths    - (2,) total widths  from and to edges of rectangle
    """
    facecolor = color

    w = [widths[1]+additional_w, widths[2]+additional_w, widths[3]+additional_w]
    # bottom left
    b = (center[1] - w[1]/2., center[2] - w[2]/2., center[3] - w[3]/2.)

    # Plot the top and bottom of the thing
    (i,j,k) = (1,2,3)
        x = [b[i], b[i], b[i]+w[i], b[i]+w[i],b[i]]
        y = [b[j], b[j]+w[j], b[j]+w[j], b[j],b[j]]
        z = b[k]*ones(5)

        fig = plot!(fig,x,y,z,linetype=:surface,colorbar=false,
            seriestype = [:shape,],
            c  = color,alpha = 0.2, seriesalpha=0.2,)

        x = [b[i], b[i], b[i]+w[i], b[i]+w[i],b[i]]
        y = [b[j], b[j]+w[j], b[j]+w[j], b[j],b[j]]

    # Why do we have this +- eps?
    # Because it will crash if all the x or y points are the same
    # which is what happens when you plot the vertical faces
    # No, i don't know why. Uncomment and see for yourself.
    eps = 1e-4
    z = b[k]*ones(5) .+ w[k]; z[1] += eps; z[2] -= eps

    fig = plot!(fig,x,y,z,linetype=:surface,colorbar=false,
        seriestype = [:shape,],
        c  = color,alpha = 0.2, seriesalpha=0.2,)


    # Plot the sides of the thing
        x = [b[1], b[1], b[1]+w[1], b[1]+w[1], b[1]]
        z = [b[3], b[3]+w[3], b[3]+w[3], b[3], b[3]]
        
    # Why do we have this +- eps?
        # Because it will crash if all the x or y points are the same
        # which is what happens when you plot the vertical faces
        # No, i don't know why. Uncomment and see for yourself.
        eps = 1e-4
        y = b[2]*ones(5); y[1] += eps; y[2] -= eps

    # Do YOU want to know why you need this twice, once with y-vector reversed?
    # It's something about the eps. It makes only a triangle instead of a square.
    # Of course, if you don't have the eps it crashes.
    # "Don't try to discover the truth."
        fig = plot!(fig,x,y,z,linetype=:surface,colorbar=false,
            seriestype = [:shape,],
            c  = color,alpha = 0.2, seriesalpha=0.2,)
        fig = plot!(fig,x,reverse(y),z,linetype=:surface,colorbar=false,
            seriestype = [:shape,],
            c  = color,alpha = 0.2, seriesalpha=0.2,)

        y = [b[2], b[2], b[2]+w[2], b[2]+w[2],b[2]]
        z = [b[3], b[3]+w[3], b[3]+w[3], b[3],b[3]]

        x = b[1]*ones(5) .+ w[1]; x[1] -= eps; x[2] += eps

    fig = plot!(fig,x,y,z,linetype=:surface,colorbar=false,
        seriestype = [:shape,],
        c  = color,alpha = 0.2, seriesalpha=0.2,)
    fig = plot!(fig,reverse(x),y,z,linetype=:surface,colorbar=false,
        seriestype = [:shape,],
        c  = color,alpha = 0.2, seriesalpha=0.2,)

    return fig
end


# From Holly's Matlab code
function plot_torus(fig; npossible_goals = 30)
    # npossible_goals: n^2 = number of points that define the torus shape
    limits = hcat(range(0,2*π, length=npossible_goals)...)
    x = Array{Float64}(undef, 3,npossible_goals*npossible_goals);

    # Generate Toroidal taskspace
    for j = 1 : npossible_goals
        for k = 1 : npossible_goals
            q = [limits[j]; limits[k]];      
            x[:,(j-1)*npossible_goals + k] = taskspace(1,3, q);# HACK l1 = 3 l2 = 1
        end
    end
    # Plot
    fig = scatter!(fig, x[1,:], x[2,:], x[3,:], primary=false,
                   markercolor=:red, markerstrokecolor=:red, markersize=2)
    return fig
end

