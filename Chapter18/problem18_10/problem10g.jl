# using Plots
# gr()

function ψ(r,a)
    return exp(-a*r)
end

function ψ2(r,a)
    return ψ(r,a)^2
end

function dψ(r,a)
    return -a*ψ(r,a)
end

function d2ψ(r,a)
    return a^2 *ψ(r,a)
end

# 湯川ポテンシャル
function V(r,α)
    return exp(-α*r)/r
end

function Hψr2(r,a,α)
    return a*r*ψ(r,a) -a^2 *r^2 *ψ(r,a)/2 + V(r,α)*ψ(r,a)*r^2
end

# [-δ,δ]の一様乱数
function point(δ)
    return δ*(2rand()-1)
end

function Metoropolis_MonteCalro(N,a,α)

    δ = 3.0
    r = 1.0

    count = 0

    # 確率分布に収束させるために最初は飛ばす
    while count < 1000

        r_trial = r + point(δ)

        w = ψ2(r_trial,a)/ψ2(r,a)

        if rand() <= w
            r = r_trial
        end

        count += 1
    end


    count = 0
    sumE = 0.0
    sumE2 = 0.0

    while count < N

        r_trial = r + point(δ)

        w = ψ2(r_trial,a)/ψ2(r,a)

        if rand() <= w
            r = r_trial
        end

        count += 1

        tmp = Hψr2(r,a,α)/(ψ(r,a)*r^2)
        sumE += tmp
        sumE2 += tmp^2

    end

    # @show count
    # @show N
    E = sumE/N
    V = sumE2/N - E^2

    println("a = $a")
    println("E = $E")
    println("V = $V")
    println()

    return V, E
end


function main()
    N = 100000

    V = 1.0
    E = 0.0

    println("start")
    α = 1.0
    for a in 1:0.1:10
        tempV, tempE = Metoropolis_MonteCalro(N,-a,α)
    # Metoropolis_MoteCalro(N,-0.14221,α)
    end

    println("finish")
end

main()
