function Integrale1(x, c, ρ, phi, psi)
    r=-phi/psi
    p=1/psi
    if psi<0
        t1=(r^2)/(p-x)
        t2=(x-p)*(c-1)^2
        t3=2*(c-1)*r*log(abs(x-p))
    else
        t1=(r^2)/(p-x)
        t2=(x-p)*(c-1)^2
        t3=2*(c-1)*r*log(abs(x-p))
    end
    return t1+t2+t3
end

function Integrale2(x, c, ρ , phi, psi)
    r=-phi/psi
    p=1/psi
    if psi<0
        t1=x*c^2
        t2=4*c*sqrt(ρ*x)
        t3=(r^2)/(x-p)
        t4=-4*r*sqrt(ρ)/(sqrt(abs(p)))*atan(sqrt(abs(x/p)))
        t5=ρ*log(x)
        t6=2*c*r*log(abs(x-p))
    else
        t1=x*c^2
        t2=4*c*sqrt(ρ*x)
        t3=(r^2)/(x-p)
        t4=4*r*sqrt(ρ/p)*atanh(sqrt(x/p))
        t5=ρ*log(x)
        t6=2*c*r*log(abs(x-p))
    end
    return t1-t2-t3+t4+t5+t6
end

function optim_c(a, b, ρ, phi, psi)
    r=-phi/psi
    p=1/psi
    t1=2*(ρ-p)
    t2=4*sqrt(ρ*a)
    t3=4*sqrt(ρ*b)
    t4=2*r*log((a-p)/(b-p))
    t6=2*(b-a-p)
    return (t1-t2+t3+t4)/t6
end
