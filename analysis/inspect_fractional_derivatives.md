```julia
include("../src/fractional_derivative.jl")
using .fractional_derivative, Plots, SpecialFunctions, LaTeXStrings, MittagLeffler
```

# Benchmark cases for the Fractional derivative

**Definitions used:**

* Grünwald–Letnikov derivative (implemente in [src/fractional_derivative.jl](../src/fractional_derivative.jl)):
    $$^{GL}D_{t}^{\alpha} = \lim_{h\rightarrow \infty} \frac{1}{h^{\alpha}} \sum_{j=0}^n (-1)^j {\alpha \choose j} f(t-jh) $$

*  Riemann-Liouville derivative:
    $$^{RL}_{a}D^{\alpha}_t f(t) = \frac{1}{\Gamma(n-\alpha)}\frac{d^n}{dt^n}\int_a^t \frac{f(\tau)}{(t-\tau)^{\alpha+1-n}}d\tau$$

Where both definitions satisfy: $^{GL}D_{t}^{\alpha} = ^{RL}_{a}D^{\alpha}_t f(t)$

**Benchmark functions:**

1. $f(t) = t^{\sigma}$; $^{RL}_{0}D^{\alpha}_t t^{\sigma} = \frac{\Gamma(\sigma + 1)}{\Gamma(\sigma - \alpha + 1)}t^{\sigma - \alpha}$

2. $f(t) = e^{t}$; $^{RL}_{0}D^{\alpha}_t e^{t} = (t- \alpha)E_{1,1-\alpha}(t-\alpha)$

Where $E_{p,q}$ is the Mittag-Leffler function: $E_{p, q} = \sum_{k=0}^{\infty} \frac{z^k}{\Gamma(p k + q)}$



```julia
function analytical1(t, sigma, alpha)
    return gamma(sigma+1)/gamma(sigma-alpha+1)*t.^(sigma - alpha)
end

function analytical2(time, alpha)
    anl = []
    for t in time
        append!(anl, [mittleff(1, 1-alpha, t)*(t^(-alpha))])
    end
    return anl    
end;
```


```julia
dt = 0.01
tf = 1
N = Int(1/dt)
times = [i*dt for i in 0:N];
```


```julia
alpha = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
sigma = 3
test_func = [(i*dt)^sigma for i in 0:N];

for (i, a) in enumerate(alpha)
    num = grunwald_letnikov(a, N, dt, test_func)
    if (i== 1)
        plot(times, analytical1(times, sigma, a), label = L"\alpha = %$a")
    else
        plot!(times, analytical1(times, sigma, a), label = L"\alpha = %$a")
    end
    if (i== length(alpha))
        plot!(times, num, color = "navy", linestyle = :dash, label = "Numerical")
    else
        plot!(times, num, color = "navy", linestyle = :dash, label = "")
    end
end
p1 = plot!(title=L"D^{\alpha}[t^{%$sigma}]");
```


```julia
alpha = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
test_func = [exp(i*dt) for i in 0:N]

for (i, a) in enumerate(alpha)
    num = grunwald_letnikov(a, N, dt, test_func)
    if (i== 1)
        plot(times, analytical2(times, a), label = L"\alpha = %$a")
    else
        plot!(times, analytical2(times, a), label = L"\alpha = %$a")
    end
    if (i== length(alpha))
        plot!(times, num, color = "navy", linestyle = :dash, label = "Numerical")
    else
        plot!(times, num, color = "navy", linestyle = :dash, label = "")
    end
end
p2 = plot!(title=L"D^{\alpha}[exp(t)]", yaxis=:log);
```


```julia
display(plot(p1, p2, size = (1000, 500)))
```

<img title="" alt="" src="inspect_fractional_derivatives_files/inspect_fractional_derivatives_6_0.png">