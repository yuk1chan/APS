# using Plots
# gr()

function ψ(x, λ)
    return exp(-λ * x^2)
end

function ψ2(x,λ)
    return exp(-2*λ * x^2)
end

function V(x)
    return (x^2)/2
end

function Hψ(x,λ)
    return -(2λ*(2(x^2)λ - 1))ψ(x,λ)/2 + V(x)ψ(x,λ)
end

# [-δ,δ]の一様乱数
function point(δ)
    return δ*(2rand()-1)
end

function Metoropolis_MonteCalro(N,λ)

    δ = 3.0
    x = 1.0

    count = 0

    # 確率分布に収束させるために最初は飛ばす
    while count < 1000

        x_trial = x + point(δ)

        w = ψ2(x_trial,λ)/ψ2(x,λ)

        if rand() <= w
            x = x_trial
            count += 1
        end
    end


    count = 0
    sumE = 0.0
    sumE2 = 0.0

    while count < N

        x_trial = x + point(δ)

        w = ψ2(x_trial,λ)/ψ2(x,λ)

        if rand() <= w
            x = x_trial
        end

        count += 1

        tmp = Hψ(x,λ)/ψ(x,λ)
        sumE += tmp
        sumE2 += tmp^2

    end

    # @show count
    # @show N
    E = sumE/N
    V = sumE2/N - E^2

    println("λ = $λ")
    println("E = $E")
    println("V = $V")
    println()
end


function main()
    N = 100000


    for λ in 0.1:0.1:1
        Metoropolis_MonteCalro(N,λ)
    end

end

main()
