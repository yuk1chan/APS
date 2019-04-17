# using Plots
# gr()

function ψ(r,a)
    return exp(-a*r)/r
end

function ψ2(r,a)
    return ψ(r,a)^2
end

function dψ(r,a)
    return -ψ(r,a)*(a+1/r)
end

function d2ψ(r,a)
    return -(a+1/r)*dψ(r,a) + ψ(r,a)/r^2
end

# 湯川ポテンシャル
function V(r,α)
    return exp(-α*r)/r
end

function Hψr2(r,a,α)
    return -(2r*dψ(r,a) + (r^2)*d2ψ(r,a))/2 + V(r,α)*ψ(r,a)*r^2
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
    if -0.0102 <= E <= -0.0101 && V > 0
        println("a = $a")
        println("E = $E")
        println("V = $V")
        println()
    end

    return V, E
end


function main()
    N = 1000000

    V = 1.0
    E = 0.0

    println("start")
    α = 1.0
    for a in 0.13:0.000001:0.15
        tempV, tempE = Metoropolis_MonteCalro(N,-a,α)

        # 収束してそうなものから
        if -0.0102 <= tempE <= -0.0101 && tempV > 0

            # 分散が最小のものを選ぶ
            if tempV < V
                E = tempE
                V = tempV
            end

        end
    end
    # Metoropolis_MoteCalro(N,-0.14221,α)

    println("finish")
end

main()

# 試行関数の形がわからない
# 確率分布に収束しているかの確認
# a = -0.142 ~ -0.143 の間くらい？？？
