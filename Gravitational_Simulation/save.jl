using LinearAlgebra
using Statistics
using JSON
using HDF5

#=
length:km, mass:kg, time:s, angle:rad(deg in json)
(interval must be a multiple of dt, tEnd must be a multiple of interval)
=#

G = 6.67430e-20 #gravitational constant (G)
AU = 1.495978e8 #AU

###### Sub Simulation Function ######
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
function calcPosition(data, t, N)
    pos = zeros(N,3)
    for i in 2:N
        e = data[i,1] #eccentricity
        IN = data[i,2] #orbital inclination
        OM = data[i,3] #ascending necliptic longitude
        W = data[i,4] #perihelion argument
        m0 = data[i,5] #Meta-period average near-point angle
        a = data[i,6] #semi-major axis
        p = data[i,7] / 86400 #Orbital period (days)

        E = calcE(calcM(m0,p,t), e) #eccentric anomaly

        #Find the coordinates (X,Y) on the orbital plane
        X = a * (cos(E) -e) #Make the origin the sun
        Y = a * sqrt(1.0 - e^2) * sin(E) 

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

        pos[i,:] = A * [
            X
            Y
            0.0
        ]
    end
    return pos
end

####### Simulation Function #######
function nBodySave(N, mass, pos, vel, tEnd, dt, chunk_size, foldername, filename)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    acc = getAcc(pos, mass, G)
    Nt = Int(ceil(tEnd/dt)) #number of steps

    f = h5open("./results/$(foldername)/$(filename).h5", "w")

    print("Simulated in")
    poses_nbody_save = zeros((N,3,chunk_size)) #Initial positions are not included
    @time for chunk in 1:Int(Nt/chunk_size)
        for i in 1:chunk_size
            vel += acc * dt/2.0
            pos += vel * dt
            acc = getAcc( pos, mass, G)
            vel += acc * dt/2.0
            poses_nbody_save[:,:,i] = pos
        end
        write(f, "/data/$chunk", poses_nbody_save)
    end

    write(f, "/params/n", N)
    write(f, "/params/tend", tEnd)
    write(f, "/params/dt", dt)
    write(f, "/params/chunksize", chunk_size)

    close(f)
end
function nBodySaveDistance(N, mass, pos, vel, tEnd, dt, chunk_size, target_index, foldername, filename)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    acc = getAcc(pos, mass, G)
    Nt = Int(ceil(tEnd/dt)) #number of steps

    f = h5open("./results/$(foldername)/$(filename).h5", "w")

    print("Simulated in")
    distance_sun_target = zeros(chunk_size)
    @time for chunk in 1:Int(Nt/chunk_size)
        for i in 1:chunk_size
            vel += acc * dt/2.0
            pos += vel * dt
            acc = getAcc( pos, mass, G)
            vel += acc * dt/2.0

            moon = 9
            # For Earth and Moon
            barycenter = [
                pos[4,1] * mass[4] + pos[moon,1] * mass[moon]
                pos[4,2] * mass[4] + pos[moon,2] * mass[moon]
                pos[4,3] * mass[4] + pos[moon,3] * mass[moon]
            ] ./ (mass[4] + mass[moon])

            distance_sun_target[i] = hypot(
                pos[1,1] - barycenter[1],
                pos[1,2] - barycenter[2],
                pos[1,3] - barycenter[3],
            )
        end
        write(f, "/data/$chunk", distance_sun_target)
    end

    write(f, "/params/n", N)
    write(f, "/params/tend", tEnd)
    write(f, "/params/dt", dt)
    write(f, "/params/chunksize", chunk_size)

    close(f)
end
function keplerSave(N, data, tEnd, dt, chunk_size, foldername, filename)
    Nt = Int(ceil(tEnd/dt)) #number of steps
    daystep = dt/86400

    f = h5open("./results/$(foldername)/$(filename).h5", "w")

    print("Simulated in")
    poses_kepler_save = zeros((N,3,chunk_size)) #Initial positions are not included
    @time for chunk in 1:Int(Nt/chunk_size)
        for i in 1:chunk_size
            poses_kepler_save[:,:,i] = calcPosition(
                data, 
                EPOCH+(chunk_size*(chunk-1)+i)*daystep,
                N
            )
        end
        write(f, "/data/$chunk", poses_kepler_save)
    end

    write(f, "/params/n", N)
    write(f, "/params/tend", tEnd)
    write(f, "/params/dt", dt)
    write(f, "/params/chunksize", chunk_size)

    close(f)
end

#Abbreviations list for the function below
abbreviations = Dict(
    "sun" => "sun", "mercury" => "mer", "venus" => "ven",
    "earth" => "ear", "mars" => "mar", "jupiter" => "jup",
    "saturn" => "sat", "uranus" => "ura", "neptune" => "nep"
)
#Make the name of result file.
function generate_filename_abbr(included_planets::Vector{String})
    abbreviated = [abbreviations[planet] for planet in included_planets]
    joined_planets = join(abbreviated, "_")
    return "simulation_$(joined_planets)"
end
#Make the name of result folder
function generate_folder_name(dt::Float64, tEnd::Float64)
    dt_str = dt >= 3600 ? "$(Int(dt / 3600))hr" : dt >= 60 ? "$(Int(dt / 60))min" : "$(Int(dt))s"
    tEnd_str = tEnd >= 86400 ? "$(Int(tEnd / 86400))d" : tEnd >= 3600 ? "$(Int(tEnd / 3600))hr" : "$(Int(tEnd))s"
    return "dt$(dt_str)_tEnd$(tEnd_str)"
end

#Load pos and vel from json
function loadVectors(key, epoch)
    file = JSON.parsefile("../Resources/nasa_data.json")
    data = file[epoch][key]
    pos = [
        data["x"],
        data["y"],
        data["z"]
    ]
    vel = [
        data["vx"],
        data["vy"],
        data["vz"]
    ]
    return pos, vel, data["mass"]
end
function loadElements(key, epoch)
    file = JSON.parsefile("../Resources/nasa_data.json")
    data = file[epoch][key]
    return [
        data["ec"],
        deg2rad(data["in"]),
        deg2rad(data["om"]),
        deg2rad(data["w"]),
        deg2rad(data["ma"]),
        data["a"],
        data["pr"]
    ]
end

function main(bodies, exclude=false)
    epoch = "2024-04-01"

    N = length(bodies) #The number of bodies
    if exclude
        N = NUMBER_OF_BODIES - N
    end
    mass = zeros(N) #List of mass
    pos = zeros((N,3)) #List of initial positions
    vel = zeros((N,3)) #List of initial velocities

    elements = zeros(N,7)
    
    # 2.Calculate initial positions and velocities of the planets
    next = 1
    included_planets = Vector{String}(undef,N)
    for i in 1:NUMBER_OF_BODIES
        is_included = PLANETS[i] in bodies
        if (is_included&!exclude)|(!is_included&exclude) # XOR
            pos[next,:], vel[next,:], mass[next] = loadVectors(PLANETS[i], epoch) # position, velocity, mass
            #elements[next,:] = loadElements(PLANETS[i], epoch) # orbital elements
            included_planets[next] = PLANETS[i]
            next += 1
        end
    end

    # 3.Simulate the motions of planets
    tEnd = 60.0 * 60.0 * 24.0 * 360 * 1000 #Endtime
    dt = 60.0 * 60.0 #Delta time
    #nBodySave(N, mass, pos, vel, tEnd, dt, 8640, generate_folder_name(dt, tEnd), "n_all_ex_moon")
    nBodySaveDistance(N, mass, pos, vel, tEnd, dt, 8640, 4 ,generate_folder_name(dt, tEnd), "n_distance_sun_earth_bary_ex_jupiter")
    #keplerSave(N, elements, tEnd, dt, 8640, generate_folder_name(dt, tEnd), "k_all")
end

EPOCH = 2460401.5
NUMBER_OF_BODIES = 10
PLANETS = ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "moon"]

if abspath(PROGRAM_FILE) == @__FILE__
    # Kepler does not support Moon
    main(
        ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "moon"],
        false
    )
end
