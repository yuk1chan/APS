
using Plots

Plots.gr()

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
   
    if x > a
        V = V0
    else
        V = 0
    end
    
    return V
end

function main()
        
    x0    = -15            # 波束の初期位置
    width = 1              # x 空間での波束の幅
    k0    = 2              # 波束の群速度
    xmax  = 20
    xmin  = -xmax
    V0    = 1
    a     = 1              # 井戸型ポテンシャルの幅の半分
    dx    = 0.4            # 刻み幅
    n     = Int64((xmax-xmin)/dx)
    dt    = 0.1            # 時間間隔
    t = 0
    tmax  = 20
    
    puls = 10
    N = n+puls
    Re = zeros(Float64,N)
    Im = zeros(Float64,N)
    Imold = zeros(Float64,N)
    V = zeros(Float64,N)
    
    initial_packet(Re,Im,x0,k0,width,xmin,n,dx,dt)
    
    
    for i in 1:n
        V[i] = calcV(i,V0,a)
    end
    
    anim = @animate for i in 1:Int64(tmax/dt)
        
        evolve(Re,Im,Imold,V0,a,dx,dt,xmin,n)
        plot(V,label="potential")
        plot!(Re,label="Real(psi)")
        plot!(Im,label="Imag(psi)")
        plot!(ylim=(-V0-0.1,V0+0.1))
        
    end
    
    gif(anim, "./anim.gif", fps = 15)

end
main()
