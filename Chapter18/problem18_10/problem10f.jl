# using Plots
# gr()

function ψ(r,a)
    return exp(-r/a)
end

function ψ2(r,a)
    return exp(-2r/a)
end

function V(r)
    return -1/r
end

function Hψr2(r,a)
    return (-1/2)*(-2ψ(r,a)/a + (r^2)ψ(r,a)/a^2) + V(r)ψ(r,a)*r^2
    # return (-2ψ(r,a)/a + (r^2)ψ(r,a)/a^2) + V(r)ψ(r,a)*r^2
end

# [-δ,δ]の一様乱数
function point(δ)
    return δ*(2rand()-1)
end

function Metoropolis_MonteCalro(N,a)

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

        tmp = Hψr2(r,a)/(ψ(r,a)*r^2)
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
end


function main()
    N = 500000

    for a in 0.1:0.1:2
        Metoropolis_MonteCalro(N,a)
    end
    # Metoropolis_MoteCalro(N,1.0)

    # a = 1.0 が理論値
end

main()
