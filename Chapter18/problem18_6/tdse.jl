
using Plots

# @animate を高速にするためのもの
Plots.gr()
ENV["PLOTS_TEST"] = "true"

function evolve(Re,Im,Imold,V0,a,dx,dt,xmin,n)
    
    for i in 2:n
        x = xmin +(i-2)*dx
    
        HIm = calcV(x,V0,a)*Im[i] - 0.5(Im[i+1] -2Im[i] +Im[i-1])/dx^2
        
        # dtの整数倍の時刻で定義された実部
        Re[i] = Re[i] + HIm*dt
    end

    for i in 2:n
        x = xmin + (i-2)*dx
        
        # 実部よりdt/2だけ前の値
        Imold[i] = Im[i]
        HRe = calcV(x,V0,a)*Re[i] -0.5(Re[i+1]-2Re[i]+Re[i-1])/dx^2
        Im[i] = Im[i] - HRe*dt
    end
    
    return Re,Im,Imold
end

function initial_packet(Re,Im,x0,k0,width,xmin,n,dx,dt)
    
    # 初期状態としてのガウス型の波束
    
    delta2 = width^2
    A = (2π*delta2)^(-0.25)
    b = k0*dt/2
    
    for i in 1:n
        x = xmin + (i-1)*dx
        arg  = 0.25*(x-x0)^2/delta2
        arg2 = 0.25*(x-x0-0.5*k0*dt)^2/delta2

        Re[i] = A*cos(k0*(x-x0))*exp(-arg)
        Im[i] = A*sin(k0*(x-x0-0.5*b))*exp(-arg2)

    end
    
    return Re,Im
end

function calcV(x,V0,a)
   
    if x >= a
        V = V0
    else
        V = 0
    end
    
    return V
end

function main()
        
    unit = 1
    x0    = -10            # 波束の初期位置
    width = 1              # x 空間での波束の幅
    k0    = 2              # 波束の群速度
    xmax  = 20
    xmin  = -xmax
    V0    = 2*unit
    a     = 1              # 井戸型ポテンシャルの幅の半分
    dx    = 0.4/sqrt(unit)            # 刻み幅
#     n     = Int64((xmax-xmin)/dx)
    n = 200
    dt    = 0.1/unit            # 時間間隔
    t = 0
    tmax  = 20
    
    plus= 10
    N = n+plus
    x = zeros(Float64,N)
    Re = zeros(Float64,N)
    Im = zeros(Float64,N)
    Imold = zeros(Float64,N)
    V = zeros(Float64,N)
    
    initial_packet(Re,Im,x0,k0,width,xmin,n,dx,dt)
    
    
    for i in 2:n
        x[i] = xmin +(i-2)*dx
        V[i] = calcV(x[i],V0,a)
    end

    unit = 10^2
    anim = @animate for i in 1:Int64(1000)
        
        evolve(Re,Im,Imold,V0,a,dx,dt,xmin,n)
        t += dt
        t = round(t*unit)/unit
        P = Re.^2 + Im.*Imold
        plot(V,label="potential")
        plot!(Re,label="Re(psi)")
        plot!(Im,label="Im(psi)")
        plot!(P ,label="|psi|")
        plot!(ylim=(-1,V0+0.1),title="t=$t")
        
        P = sum(Re.*Re + Im.*Imold)
        P = P*dx
        ψ = Re.*Re + Im.*Imold
        ψ = ψ./P
        
#         println("P = $P")
#         println(sum(ψ)*dx)
        
        temp1 = sum(x.^2 .*ψ)*dx
        temp2 = (sum(x.*ψ)^2)*dx*dx
        temp3 = sum(x.*ψ)*dx
        w = temp1 - temp2 # 確率密度分布の分散
#         println("temp1 = $temp1 \t temp2 = $temp2")
        println("xの平均値: $temp3 \t t = $t")
    end
    
    w = width
    name = "k0=$k0,w=$w"
    
    gif(anim, "./$name.gif", fps = 60)
    

end
main()
