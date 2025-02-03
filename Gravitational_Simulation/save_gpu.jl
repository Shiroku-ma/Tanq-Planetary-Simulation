using LinearAlgebra
using Statistics
using JSON
using HDF5
using Metal

include("save.jl")

function mat_trans_dif(mt_pos, diff)
    i, j = thread_position_in_grid_2d()  # 各スレッドの 2D 座標を取得
    @inbounds begin
        diff[i,j,1] = mt_pos[i,1] - mt_pos[j,1]
        diff[i,j,2] = mt_pos[i,2] - mt_pos[j,2]
        diff[i,j,3] = mt_pos[i,3] - mt_pos[j,3]
    end
    return
end

function mat_univ_grav(mt_mass, diff, acc)
    i, j = thread_position_in_grid_2d()  # 各スレッドの 2D 座標を取得
    @inbounds begin
        inv_r3 = (diff[i,j,1]^2 + diff[i,j,2]^2 + diff[i,j,3]^2 + Float32(1.0e-15))^Float32(-1.5)
        acc[i,j,1] = diff[i,j,1] * inv_r3 * mt_mass[i] * Float32(6.67430e-20)
        acc[i,j,2] = diff[i,j,2] * inv_r3 * mt_mass[i] * Float32(6.67430e-20)
        acc[i,j,3] = diff[i,j,3] * inv_r3 * mt_mass[i] * Float32(6.67430e-20)
    end
    return
end

function copy_to_buff(s, poses_nbody_save, pos)
    i, j = thread_position_in_grid_2d() 
    poses_nbody_save[i,j,s] = pos[i,j]
    return
end

function scale_and_add(a, mat1, mat2)
    i, j = thread_position_in_grid_2d()
    mat1[i,j] += mat2[i,j] * a
    return
end

function getAcc(N, mt_pos, mt_mass , G)
    dims=(N,N)
    diff = Metal.zeros(N,N,3)
    @time Metal.@sync @metal threads=dims mat_trans_dif(mt_pos, diff)

    acc_div = Metal.zeros(N,N,3)
    @time Metal.@sync @metal threads=dims mat_univ_grav(mt_mass, diff, acc_div)
    mt_acc = MtlArray(sum(acc_div, dims=1))

    #=
    For CPU
    #Distance between celestial bodies
    dx = zeros(N,N)
    dy = zeros(N,N)
    dz = zeros(N,N)
    @inbounds for j in 1:N
        @simd for i in 1:N
            dx[i,j] = pos[i,1] - pos[j,1]
            dy[i,j] = pos[i,2] - pos[j,2]
            dz[i,j] = pos[i,3] - pos[j,3]
        end
    end

    #Inverse of the cube of the distance
    ax = zeros(N)
    ay = zeros(N)
    az = zeros(N)
    @inbounds for j in 1:N
        @simd for i in 1:N
            inv_r3_gm = (dx[i,j]^2 + dy[i,j]^2 + dz[i,j]^2 + 1.0e-15)^(-1.5) * G * mass[i]
            ax[j] += dx[i,j] * inv_r3_gm
            ay[j] += dy[i,j] * inv_r3_gm
            az[j] += dz[i,j] * inv_r3_gm
        end
    end
    =#

    return mt_acc
end

####### Simulation Function #######
function nBodySave(N, mass, pos, vel, tEnd, dt, chunk_size, foldername, filename)
    Nt = Int(ceil(tEnd/dt)) #number of steps

    mt_mass = MtlArray(mass)
    mt_pos = MtlArray(pos)
    vel .-= mean(mass .* vel) / mean(mass) #center of mass system
    mt_vel = MtlArray(vel)
    mt_acc = Metal.zeros(N,3)
    mt_acc = getAcc(N, mt_pos, mt_mass, G)

    #=
    f = h5open("./results/$(foldername)/$(filename).h5", "w")

    print("Simulated in")
    poses_nbody_save = Metal.zeros(N,3,chunk_size)
    @time for chunk in 1:Int(Nt/chunk_size)
        for i in 1:chunk_size
            Metal.@sync @metal threads=(N,3) scale_and_add(dt/2.0, mt_vel, mt_acc)
            Metal.@sync @metal threads=(N,3) scale_and_add(dt, mt_pos, mt_vel)
            mt_acc = getAcc(N,mt_pos, mt_mass, G)
            Metal.@sync @metal threads=(N,3) scale_and_add(dt/2.0, mt_vel, mt_acc)
            Metal.@sync @metal threads=(N,3) copy_to_buff(i, poses_nbody_save, mt_pos)
        end
        write(f, "/data/$chunk", unsafe_wrap(Array, poses_nbody_save, (N,3,chunk_size)))
    end

    write(f, "/params/n", N)
    write(f, "/params/tend", tEnd)
    write(f, "/params/dt", dt)
    write(f, "/params/chunksize", chunk_size)

    close(f)
    =#
end

function main(bodies, exclude=false)
    epoch = "2024-04-01"

    N = length(bodies) #The number of bodies
    if exclude
        N = NUMBER_OF_BODIES - N
    end
    mass = zeros(Float32, N) #List of mass
    pos = zeros(Float32, (N,3)) #List of initial positions
    vel = zeros(Float32, (N,3)) #List of initial velocities
    
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
    tEnd = 60.0 * 60.0 * 24.0 * 360 * 10 #Endtime
    dt = 60.0 * 60.0 #Delta time
    nBodySave(N, mass, pos, vel, tEnd, dt, 8640, "tmp", "gputest")
end

G = Float32(6.67430e-20)
EPSILON = Float32(1.0e-15)
INV3 = Float32(-1.5)

main(
    ["sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "moon"],
    false
)
