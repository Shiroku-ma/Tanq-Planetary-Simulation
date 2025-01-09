using LinearAlgebra
using Statistics
using JSON
using HDF5

#=
length:km, mass:kg, time:s, angle:rad(deg in json)
(interval must be a multiple of dt, tEnd must be a multiple of interval)
=#

colors = [:orange, :black, :gold, :blue, :red, :brown, :bisque4, :gray, :blue] #惑星の色
G = 6.67430e-20 #gravitational constant (G)
AU = 1.495978e8 #AU

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
    # 時間間隔(dt)のフォーマット（秒を分・時間に変換）
    dt_str = dt >= 3600 ? "$(Int(dt / 3600))hr" : dt >= 60 ? "$(Int(dt / 60))min" : "$(Int(dt))s"
    # 総時間(tEnd)のフォーマット（秒を日・時間に変換）
    tEnd_str = tEnd >= 86400 ? "$(Int(tEnd / 86400))d" : tEnd >= 3600 ? "$(Int(tEnd / 3600))hr" : "$(Int(tEnd))s"
    return "dt$(dt_str)_tEnd$(tEnd_str)"
end

#Load pos and vel from json
function loadData(key, epoch)
    file = JSON.parsefile("../Resources/nasa_pos_data.json")
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

function main(bodies, exclude=false)
    EPOCH = "2024-04-01"

    N = length(bodies) #The number of bodies
    if exclude
        N = NUMBER_OF_BODIES - N
    end
    mass = zeros(N) #List of mass
    pos = zeros((N,3)) #List of initial positions
    vel = zeros((N,3)) #List of initial velocities
    
    # 2.Calculate initial positions and velocities of the planets
    next = 1
    included_planets = Vector{String}(undef,N)
    for i in 1:NUMBER_OF_BODIES
        is_included = PLANETS[i] in bodies
        if (is_included*!exclude+!is_included*exclude) == 1 # XOR
            pos[next,:], vel[next,:], mass[next] = loadData(PLANETS[i], EPOCH)
            included_planets[next] = PLANETS[i]
            next += 1
        end
    end

    # 3.Simulate the motions of planets
    tEnd = 60.0 * 60.0 * 24.0 * 360 * 10 #Endtime
    dt = 60.0 * 60.0 #Delta time
    nBodySave(N, mass, pos, vel, tEnd, dt, 8640, generate_folder_name(dt, tEnd), generate_filename_abbr(included_planets))
end

NUMBER_OF_BODIES = 9
PLANETS = ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"]

main(
    ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"],
    false
)
