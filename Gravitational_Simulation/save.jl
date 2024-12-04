using LinearAlgebra
using Statistics
using JSON
using HDF5

#=
length:km, mass:kg, time:s, angle:rad(deg in json)
(interval must be a multiple of dt, tEnd must be a multiple of interval)
=#

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
function nBodySave(N, mass, pos, vel, tEnd, dt, chunk_size)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    acc = getAcc(pos, mass, G)
    Nt = Int(ceil(tEnd/dt)) #number of steps

    f = h5open("./results/dt1hr_tEnd1080d/all.h5", "w")

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
#Load pos and vel from json
function loadData(key, epoch, sun)
    file = JSON.parsefile("../Resources/nasa_pos_data.json")
    if sun
        return file[epoch]["sun_m"]
    end
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

function main()
    EPOCH = "2024-04-01"

    N = 9 #The number of planets + 1 (sun)
    mass = zeros(N) #List of mass
    mass[1] = loadData("", EPOCH, true)
    pos = zeros((N,3)) #List of initial positions
    vel = zeros((N,3)) #List of initial velocities
    

    # 2.Calculate initial positions and velocities of the planets
    for i in 2:N
        pos[i,:], vel[i,:], mass[i] = loadData(NAMEINDEX[i], EPOCH, false)
    end

    # 3.Simulate the motions of planets
    tEnd = 60.0 * 60.0 * 24.0 * 360 * 3 #Endtime
    dt = 60.0 * 60.0 #Delta time
    nBodySave(N, mass, pos, vel, tEnd, dt, 8640)
end

main()
