using LinearAlgebra
using Statistics
using Plots
using FFMPEG
using JSON

#=
length:m, mass:kg, time:s, angle:rad(deg in json)
(interval must be a multiple of dt, tEnd must be a multiple of interval)
=#


global EPOCH
global SUN_M
NAMEINDEX = Dict(
    1 => "sun",
    2 => "mercury",
    3 => "venus",
    4 => "earth",
    5 => "mars",
    6 => "jupiter",
    7 => "saturn",
    8 => "uranus",
    9 => "neptune"
)
colors = [:orange, :black, :gold, :blue, :red, :brown, :bisque4, :gray, :blue] #惑星の色
G = 6.67430e-11 #gravitational constant (G)
AU = 1.495978e11 #AU

# <---Using Kepler's Law to find position and velocity--->
#Find the mean perigee angle (M)
function calcM(m0, p, t)
    M = m0 + 2π * (t - EPOCH) / p
    return M
end
#Find the eccentric proximal angle (E)
function calcE(M, e)
    E = M
    # Newton's method
    for i in 1:5
        E = E - (M - E + e * sin(E)) / (e * cos(E) - 1.0)
    end
    return E
end
# Find position and velocity
function calcPosAndVel(data, t)
    IN = data["IN"] #orbital inclination
    OM = data["OM"] #ascending necliptic longitude
    W = data["W"] #perihelion argument
    a = data["A"] #semi-major axis
    e = data["EC"] #eccentricity
    m0 = data["MA"] #Meta-period average near-point angle
    p = data["PR"] #Orbital period (days)

    E = calcE(calcM(m0,p,t), e) #eccentric anomaly

    #Find the coordinates (X,Y) on the orbital plane
    eee = 1.0 - e^2
    X = a * (cos(E) -e) #Make the origin the sun
    Y = a * sqrt(eee) * sin(E) 

    r = hypot(X,Y) #Distance from the Sun
    v = sqrt(G * SUN_M * (2/r - 1/a)) #orbital velocity

    #Decompose orbital velocity into xy components
    oX = X + a*e
    c = 1.0 / hypot(Y, eee * oX)
    VX = -v * Y * c
    VY = v * oX * eee * c

    #Convert orbital plane coordinates to heliocentric ecliptic coordinates
    A = [
        cos(OM) -sin(OM) 0.0
        sin(OM) cos(OM) 0.0
        0.0 0.0 1.0
    ] * [
        1.0 0.0 0.0
        0.0 cos(IN) sin(IN)
        0.0 -sin(IN) cos(IN)
    ] * [
        cos(W) -sin(W) 0.0
        sin(W) cos(W) 0.0
        0.0 0.0 1.0
    ]

    pos = A * [
        X
        Y
        0.0
    ]
    vel = A* [
        VX
        VY
        0.0
    ]

    return pos, vel
end

# <---N body simulation--->
# get acc
function getAcc(pos, mass , G)
    x = pos[:,1,1]
    y = pos[:,2,1]
    z = pos[:,3,1]

    #Distance between celestial bodies
    dx = x' .- x
    dy = y' .- y
    dz = z' .- z

    #Inverse of the cube of the distance
    inv_r3 = (dx.^2 .+ dy.^2 .+ dz.^2 .+ 1.0e-15).^(-1.5)

    #acc
    ax = G * (dx .* inv_r3) * mass
    ay = G * (dy .* inv_r3) * mass
    az = G * (dz .* inv_r3) * mass

    return hcat(ax,ay,az)
end
function nBody(N, mass, pos, vel, tEnd, dt)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    acc = getAcc(pos, mass, G)
    Nt = Int(ceil(tEnd/dt)) #number of steps
    
    poses_nbody_save = zeros((N,3,Nt+1))
	poses_nbody_save[:,:,1] = pos

    print("Simulated in")
    @time for i in 1:Nt
        #Leapfrog methodc
		vel += acc * dt/2.0
		pos += vel * dt
		acc = getAcc( pos, mass, G)
		vel += acc * dt/2.0

        poses_nbody_save[:,:,i+1] = pos
    end

    return poses_nbody_save

end

# Plot with lines
function plotAll(N, pos_save, tEnd, dt, range, title)
    p = scatter([], [], xlims=(-range,range), ylims=(-range,range), title=title, xlabel="x(au)", ylabel="y(au)", aspect_ratio=:equal, legend=false)
    for i in 2:N
        X = pos_save[i,1,:] ./ AU
        Y = pos_save[i,2,:] ./ AU
        plot!(p, X, Y, color=colors[i])
    end
    return p
end
# Plot with dots
function scatterAll(N, pos_save, tEnd, dt, range, interval, animate, title)
    p = scatter([], [], xlims=(-range,range), ylims=(-range,range), title=title, xlabel="x", ylabel="y", aspect_ratio=:equal, legend=false)
    step = Int(interval / dt)
    pos_save[:,:,:] ./= AU
    if animate
        anim = @animate for frame in 1:Int(tEnd / interval)
            for i in 1:N
                x = pos_save[i,1,frame*step]
                y = pos_save[i,2,frame*step]
                scatter!(p, [x], [y], color=colors[i], markersize=2, marker_stroke=nothing, markerstrokewidth=0)
            end
        end

        return p, anim
    else
        for frame in 1:Int(tEnd / interval)
            for i in 1:N
                x = pos_save[i,1,frame*step]
                y = pos_save[i,2,frame*step]
                scatter!(p, [x], [y], color=colors[i], markersize=2, marker_stroke=nothing, markerstrokewidth=0)
            end
        end
        return p, nothing
    end
end
#Load orbits data, and convert from degree to radian
function loadData(epoch)
    file = JSON.parsefile("../Resources/nasa_data.json")
    planets = file[epoch]
    for (key, _s) in planets
        if key != "epoch" && key != "G" && key != "sun_m" 
            for (e, value) in planets[key]
                if e == "IN" || e == "OM" || e == "W" || e == "MA"
                    planets[key][e] = deg2rad(value)
                end
            end
        end
    end
    return planets
end


function main()
    global EPOCH
    global SUN_M

    # 1.Load orbit element data from json file and convert degrees to radians
    planets = loadData("2024-04-01")
    EPOCH = planets["epoch"]
    SUN_M = planets["sun_m"]
    

    N = 6 #The number of planets + 1 (sun)
    mass = zeros(N) #List of mass
    mass[1] = SUN_M
    pos = zeros((N,3)) #List of initial positions
    vel = zeros((N,3)) #List of initial velocities
    

    # 2.Calculate initial positions and velocities of the planets
    for i in 2:N
        data = planets[NAMEINDEX[i]]
        mass[i] = data["mass"]
        pos[i,:], vel[i,:] = calcPosAndVel(data, EPOCH)
    end

    # 3.Simulate the motions of planets
    tEnd = 60.0 * 60.0 * 24.0 * 360 #Endtime
    dt = 60.0 * 60.0 #Delta time
    poses_nbody_save = nBody(N, mass, pos, vel, tEnd, dt)

    # 4.Plot
    p1 = plotAll(N, poses_nbody_save, tEnd, dt, 5.4, "Planet Orbits")
    scatter_interval = 864000
    p2, anim2 = scatterAll(N, poses_nbody_save, tEnd, dt, 5.4, scatter_interval, true, "Planet Motinos")
    display(p1)
    display(p2)
    gif(anim2, fps=10, "./images/$N-$tEnd-$dt-$scatter_interval.gif")
end

gr()
main()
