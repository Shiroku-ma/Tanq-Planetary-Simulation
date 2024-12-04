using LinearAlgebra
using Statistics
using Plots
using FFMPEG
using JSON
using HDF5

#=
length:km, mass:kg, time:s, angle:rad(deg in json)
(interval must be a multiple of dt, tEnd must be a multiple of interval)
=#

colors = [:orange, :black, :gold, :blue, :red, :brown, :bisque4, :gray, :blue] #惑星の色
AU = 1.495978e8 #AU

function plotAllWithChunks(N, file, tEnd, dt, range, title, chunk_size)
    p = scatter([], [], xlims=(-range,range), ylims=(-range,range), title=title, xlabel="x(au)", ylabel="y(au)", aspect_ratio=:equal, legend=false)
    Nt = Int(ceil(tEnd/dt))
    for c in 1:Int(Nt/chunk_size)
        chunk = read(file, "/data/$c")
        for i in 2:N
            X = chunk[i,1,:] ./ AU
            Y = chunk[i,2,:] ./ AU
            plot!(p, X, Y, color=colors[i])
        end
    end
    return p 
end

function main()
    f = h5open(".results/dt1hr_tEnd1080d/all.h5", "r")
    N = read(f, "/params/n")
    tEnd = read(f, "/params/tend")
    dt = read(f, "/params/dt")
    chunk_size = read(f, "/params/chunksize")

    @time p1 = plotAllWithChunks(N, f, tEnd, dt, 5.4, "Planets Orbits", chunk_size)

    close(f)
    display(p1)
end
gr()
main()
