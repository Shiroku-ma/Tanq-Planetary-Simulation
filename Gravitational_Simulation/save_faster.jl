using LinearAlgebra
using Statistics
using JSON
using HDF5

include("save.jl")

function getAcc(N, pos, mass, dpos, acc, G)
    fill!(acc, 0.0)

    #Distance between celestial bodies
    @inbounds for dim in 1:3
        @inbounds for j in 1:N-1
            @inbounds for i in j+1:N
                temp = pos[i,dim] - pos[j,dim]
                dpos[i,j,dim] = temp
                dpos[j,i,dim] = -temp
            end
        end
    end

    #Inverse of the cube of the distance
    @inbounds for j in 1:N
        @inbounds for i in 1:N
            inv_r3_gm = (dpos[i,j,1]^2 + dpos[i,j,2]^2 + dpos[i,j,3]^2 + 1.0e-15)^(-1.5) * G * mass[i]
            @inbounds for dim in 1:3
                acc[j,dim] += dpos[i,j,dim] * inv_r3_gm
            end
        end
    end

    return acc
end

####### Simulation Function #######
function nBodySave(N, mass, pos, vel, tEnd, dt, chunk_size, foldername, filename)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    acc = getAcc(N, pos, mass, G)
    Nt = Int(ceil(tEnd/dt)) #number of steps

    f = h5open("./results/$(foldername)/$(filename).h5", "w")

    poses_nbody_save = zeros((N,3,chunk_size)) #Initial positions are not included
    @time for chunk in 1:Int(Nt/chunk_size)
        for i in 1:chunk_size
            for j in eachindex(acc)
                vel[j] += acc[j] * dt/2.0
                pos[j] += vel[j] * dt
            end
            acc = getAcc(N,pos, mass, G)
            for j in eachindex(acc)
                vel[j] += acc[j] * dt/2.0
            end
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
function nBodySaveDistance(N, mass, pos, vel, tEnd, dt, chunk_size, target_index, moon_index, foldername, filename)
    dpos = zeros(N,N,3)
    acc = zeros(N,3)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    acc = getAcc(N, pos, mass, dpos, acc, G)
    Nt = Int(ceil(tEnd/dt)) #number of steps

    f = h5open("./results/$(foldername)/$(filename).h5", "w")

    distance_sun_target = zeros(chunk_size)
    barycenter = zeros(3)
    barymass = 0.0
    @time for chunk in 1:Int(Nt/chunk_size)
        for i in 1:chunk_size
            for j in eachindex(acc)
                vel[j] += acc[j] * dt/2.0
                pos[j] += vel[j] * dt
            end
            acc = getAcc(N, pos, mass, dpos, acc, G)
            for j in eachindex(acc)
                vel[j] += acc[j] * dt/2.0
            end

            # Calculate distance between the sun and the barycenter of a planet and its moon
            # The planet and its moon are designated by index
            # moon_index should be the same as target_index when the moon is excluded
            barymass = mass[target_index] + mass[moon_index]
            barycenter[1] = (pos[target_index,1] * mass[target_index] + pos[moon_index,1] * mass[moon_index]) / barymass
            barycenter[2] = (pos[target_index,2] * mass[target_index] + pos[moon_index,2] * mass[moon_index]) / barymass
            barycenter[3] = (pos[target_index,3] * mass[target_index] + pos[moon_index,3] * mass[moon_index]) / barymass
            distance_sun_target[i] = hypot(
                pos[1,1] - barycenter[1],
                pos[1,2] - barycenter[2],
                pos[1,3] - barycenter[3],
            )
        end
        write(f, "/data/" * string(chunk), distance_sun_target)
    end

    write(f, "/params/n", N)
    write(f, "/params/tend", tEnd)
    write(f, "/params/dt", dt)
    write(f, "/params/chunksize", chunk_size)

    close(f)
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
    
    # 2.Calculate initial positions and velocities of the planets
    next = 1
    included_planets = Vector{String}(undef,N)
    for i in 1:NUMBER_OF_BODIES
        is_included = PLANETS[i] in bodies
        if (is_included&!exclude)|(!is_included&exclude) # XOR
            pos[next,:], vel[next,:], mass[next] = loadVectors(PLANETS[i], epoch) # position, velocity, mass
            included_planets[next] = PLANETS[i]
            next += 1
        end
    end

    # 3.Simulate the motions of planets
    tEnd = 60.0 * 60.0 * 24.0 * 365.25 * 500000 #Endtime
    dt = 60.0 * 60.0 #Delta time

    #nBodySave(N, mass, pos, vel, tEnd, dt, 86400, generate_folder_name(dt, tEnd), "n_all")
    nBodySaveDistance(
        N, mass, pos, vel, tEnd, dt,
        876600, findfirst(x->x=="earth", included_planets), findfirst(x->x=="moon", included_planets),
        "tmp", "n_distance_sun_earth_bary"
        )
    #keplerSave(N, elements, tEnd, dt, 86400, generate_folder_name(dt, tEnd), "k_all_bary")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(
        [],
        true
    )
end
