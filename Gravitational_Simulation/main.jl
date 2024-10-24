using LinearAlgebra
using Statistics
using Plots
using FFMPEG
using JSON

#=
長さ:m, 質量:kg, 時間:s, 角度:rad(jsonはdeg)
intervalはdtの倍数、tEndはintervalの倍数でなければならない
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
colors = [:blue, :black, :gold, :blue, :red, :brown, :bisque4, :gray, :blue] #惑星の色
G = 6.67430e-11 #万有引力定数
AU = 1.495978e11 #天文単位

# <---ケプラーの法則を利用して位置と速度を求める--->
#平均近点角(M)を求める
function calcM(m0, p, t)
    M = m0 + 2π * (t - EPOCH) / p
    return M
end
#離心近点角(E)を求める
function calcE(M, e)
    E = M
    # ニュートン法
    for i in 1:5
        E = E - (M - E + e * sin(E)) / (e * cos(E) - 1.0)
    end
    return E
end
# 位置と速度を求める
function calcPosAndVel(data, t)
    IN = data["IN"] #軌道傾斜角
    OM = data["OM"] #昇交点黄経
    W = data["W"] #近日点引数
    a = data["A"] #軌道長半径
    e = data["EC"] #離心率
    m0 = data["MA"] #元期平均近点角
    p = data["PR"] #公転周期(日)

    E = calcE(calcM(m0,p,t), e) #離心近点角

    #軌道面上での座標(X,Y)を求める
    eee = 1.0 - e^2
    X = a * (cos(E) -e) #原点を太陽にする
    Y = a * sqrt(eee) * sin(E) 

    r = hypot(X,Y) #太陽との距離
    v = sqrt(G * SUN_M * (2/r - 1/a)) #軌道速度

    #軌道速度をxy成分に分解
    oX = X + a*e #原点を中心に戻す
    c = 1.0 / hypot(Y, eee * oX)
    VX = -v * Y * c
    VY = v * oX * eee * c

    #軌道面座標を日心黄道座標に変換
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

# <---N体シミュレーション--->
# 加速度を取得
function getAcc(pos, mass , G)
    x = pos[:,1,1]
    y = pos[:,2,1]
    z = pos[:,3,1]

    #天体間の距離
    dx = x' .- x
    dy = y' .- y
    dz = z' .- z

    #距離の3乗の逆数
    inv_r3 = (dx.^2 .+ dy.^2 .+ dz.^2 .+ 1.0e-15).^(-1.5)

    #加速度
    ax = G * (dx .* inv_r3) * mass
    ay = G * (dy .* inv_r3) * mass
    az = G * (dz .* inv_r3) * mass

    #連結
    return hcat(ax,ay,az)
end
function nBody(N, mass, pos, vel)
    tEnd = 60.0 * 60.0 * 24.0 * 360 #終了時間
    dt = 60.0 * 60.0
    vel .-= mean(mass .* vel) / mean(mass) #重心系
    acc = getAcc(pos, mass, G)
    Nt = Int(ceil(tEnd/dt)) #ステップ数
    
    #行列　[天体数、 3d、 Nt+1]
    poses_nbody_save = zeros((N,3,Nt+1))
	poses_nbody_save[:,:,1] = pos

    @time for i in 1:Nt
        #リープ・フロッグ法
		vel += acc * dt/2.0
		pos += vel * dt
		acc = getAcc( pos, mass, G)
		vel += acc * dt/2.0

        poses_nbody_save[:,:,i+1] = pos
    end

    @time p1 = plotAll(N, poses_nbody_save, tEnd, dt, 5.4)
    @time p2, anim2 = scatterAll(N, poses_nbody_save, tEnd, dt, 5.4, 864000, true)

    @time display(p1)
    @time display(p2)
    gif(anim2, fps=10, "images/scatter.gif")
    #savefig("images/plot.png")

end

# 全てプロット
function plotAll(N, pos_save, tEnd, dt, range)
    p = scatter([], [], xlims=(-range,range), ylims=(-range,range), title="Positions", xlabel="x(au)", ylabel="y(au)", aspect_ratio=:equal, legend=false)
    for i in 2:N
        X = pos_save[i,1,:] ./ AU
        Y = pos_save[i,2,:] ./ AU
        plot!(p, X, Y, color=colors[i])
    end
    return p
end
# 全て点でプロット
function scatterAll(N, pos_save, tEnd, dt, range, interval, animate)
    p = scatter([], [], xlims=(-range,range), ylims=(-range,range), title="Positions(All)", xlabel="x", ylabel="y", aspect_ratio=:equal, legend=false)
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
#元期軌道要素のデータを読み込み、ラジアンに変換
function loadData(epoch)
    file = JSON.parsefile("./nasa_data.json")
    planets = file[epoch]
    for (key, _s) in planets
        if key != "epoch" && key != "G" && key != "sun_m" # 惑星以外の項目は除く
            for (e, value) in planets[key] # 各惑星について
                if e == "IN" || e == "OM" || e == "W" || e == "MA" # 角度以外の項目は除く
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

    planets = loadData("2024-04-01")
    EPOCH = planets["epoch"]
    SUN_M = planets["sun_m"]

    N = 6
    mass = zeros(N)
    mass[1] = SUN_M
    pos = zeros((N,3))
    vel = zeros((N,3))

    for i in 2:N
        data = planets[NAMEINDEX[i]]
        mass[i] = data["mass"]
        pos[i,:], vel[i,:] = calcPosAndVel(data, EPOCH)
    end
    #=
    開始時刻は元期でなくてもいいが、誤差を考えると元期に近い方が望ましい
    =#

    nBody(N, mass, pos, vel)
end

gr() #Plotのバックエンドを指定（デフォルト）
main()
