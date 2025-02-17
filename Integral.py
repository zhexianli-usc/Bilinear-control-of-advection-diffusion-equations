import numpy as np

from scipy.integrate import tplquad, dblquad, quad, nquad
from scipy.special import erfc


# Define the integrand f(z, y, x).
# Here f(x,y,z) = x + y + z, so we reorder arguments as (z, y, x).
i = 1j
T = 2
L = 1
ray_right  = np.exp(i * np.pi / 8)
ray_left = np.exp(7 * i * np.pi / 8)

def Delta(k):
    return np.exp(i * k * L) - np.exp(-i * k * L)

def omega(k, v):
    return k ** 2 + i * v * k

def StateHatInitial(k):
    return -i * (1 - np.exp(-1 - i * k)) / (k - i)

def StateInitialDeri(k):
    return - np.exp(-1 - i * k) * (2 - np.exp(1 + i * k) + i * k) / (- i + k) ** 2

def Boundary_left(k, t):
    value = (np.exp(k * t) - 1) / k
    if np.isnan(value):
        return t
    else:
        return (np.exp(k * t) - 1) / k

def Boundary_derivative_left(k, t, v):
    value = i * (1 + np.exp(k * t * (k + i * v)) * (-1 + k**2 * t + i * k * v)) / (k * (k + i * v) ** 2)
    if np.isnan(value):
        return 0
    else:
        return value

def integrand_state_real(k, x, t, v):
    return np.real( np.exp(i * k * x - omega(k, v) * t) * StateHatInitial(k))

def integrand_state_complex(r, x, t, v, ray):
    k_trans = - r * ray - i * v
    k = r * ray
    
    return - i * (2 * k + i * v) / (k + i * v) * np.exp(i * k * x) * (1 - np.exp(- omega(k,v) * t)) / k
    # return np.real( - np.exp( - omega(k_ori, v) * t) * 
    #                ( 2 * np.sin(k_ori * x) * i * np.exp(i * k_ori * L) * StateHatInitial(k_ori) +
    #                  2 * np.sin(k_ori * (L - x)) * i * StateHatInitial(k_trans)
    #                 )
    #                )

def integrand_state_complex_sum(r, x,t, v):
    return np.real(integrand_state_complex(r, x, t, v, ray_right) * ray_right - integrand_state_complex(r, x, t, v, ray_left) * ray_left) / (2*np.pi)

def integrand_de_real(k, x, t, v):
    return np.real(- i * k * t * np.exp(i * k * x - omega(k, v) * t) * StateHatInitial(k))

def integrand_de_complex(r, x, t, v, ray):
    k_trans = - r * ray - i * v
    k   = r * ray
    
    exponential_term = np.exp(i * k * x - omega(k, v) * t)
    return -(np.exp(i * k * x) + exponential_term * (-1 -2 * k **2 * t - 3 * i * k * t * v + t * v**2)) / (k + i * v)**2
    
    # integrand_state  = (- np.exp( - omega(k_ori, v) * t) / Delta(k_ori) * 
    #                 ( 2 * np.sin(k_ori * x) * i * np.exp(i * k_ori * L) * StateHatInitial(k_ori) +
    #                  2 * np.sin(k_ori * (L - x)) * i * StateHatInitial(k_trans)
    #                 ))
    # #print(integrand_state)
    # return np.real( - i * k_ori * t * integrand_state - 
    #                ( np.exp( - omega(k_ori, v) * t) / Delta(k_ori) * 
    #                 (2 * np.sin(k_ori * (L - x)) * i * StateInitialDeri(k_trans) 
    #                 ) 
    #                )
    #               )

def integrand_de_complex_sum(r, x, t, v):
    return np.real(integrand_de_complex(r, x, t, v, ray_right) * ray_right - integrand_de_complex(r, x, t, v, ray_left) * ray_left) / (2 * np.pi)

def Laplace_solution(x,v,t):
    return 1/2 * np.exp(v * x /2) * ( np.exp(-v * x/2) * erfc((x - v*t) / (2 * np.sqrt(t))) + np.exp(v * x / 2) * erfc((x + v*t) / (2 * np.sqrt(t))))

# for v in range(1,10):
#     integrand = lambda r, x: integrand_state_complex_sum(r, x, 1, v)
#     print(dblquad(
#         integrand,
#         0,           # x lower limit
#         np.inf,          # x upper limit
#         lambda x: 0,      # k lower limit (can depend on x)
#         lambda x: np.inf,      # k upper limit (can depend on x)
#         #epsabs=1e-2,   # <-- Increased absolute tolerance
#         #epsrel=1e-2,    # <-- Increased relative tolerance
#                 )[0]
#         )
#     integrand = lambda x: Laplace_solution(x,v,1)
#     print(f"Laplace solution: {quad(integrand,0,100)[0]}")

# Use tplquad to integrate over x in [0,1], y in [0,2], z in [0,3].
# tplquad(func, x_min, x_max, y_min, y_max, z_min, z_max)
t_division = 10
v = np.zeros(t_division + 1) + 3
for step in range(100):
    v_grad = np.zeros(t_division + 1)
    obj = np.zeros(t_division + 1)

    for t_step in range(t_division + 1):
        # integrand = lambda k, x: integrand_de_real(k, x, T / t_division * t_step, v[t_step])
        # v_grad[t_step] += dblquad(
        #     integrand,
        #     0,           # x lower limit
        #     L,       # x upper limit
        #     lambda x: -np.inf,      # k lower limit (can depend on x)
        #     lambda x: np.inf,      # k upper limit (can depend on x)
        #     epsabs=1e-2,   # <-- Increased absolute tolerance
        #     epsrel=1e-2,    # <-- Increased relative tolerance
        # )[0]
        
        integrand = lambda r, x: integrand_de_complex_sum(r, x, T / t_division * t_step, v[t_step])
        
        
        v_grad[t_step] += \
        - dblquad(
            integrand,
            0,           # x lower limit
            100,       # x upper limit
            lambda x: 0,      # k lower limit (can depend on x)
            lambda x: np.inf,      # k upper limit (can depend on x)
            #epsabs=1e-2,   # <-- Increased absolute tolerance
            #epsrel=1e-2,    # <-- Increased relative tolerance
        )[0] \
        + v[t_step]
        
        # integrand = lambda k, x: integrand_state_real(k, x, T / t_division * t_step, v[t_step])
        
        # obj[t_step] += dblquad(
        #     integrand,
        #     0,           # x lower limit
        #     L,          # x upper limit
        #     lambda x: -np.inf,      # k lower limit (can depend on x)
        #     lambda x: np.inf,      # k upper limit (can depend on x)
        #     epsabs=1e-2,   # <-- Increased absolute tolerance
        #     epsrel=1e-2,    # <-- Increased relative tolerance
        # )[0]
        
        integrand = lambda r, x: integrand_state_complex_sum(r, x, T / t_division * t_step, v[t_step])
        
        obj[t_step] += \
        - dblquad(
            integrand,
            0,           # x lower limit
            np.inf,          # x upper limit
            lambda x: 0,      # k lower limit (can depend on x)
            lambda x: np.inf,      # k upper limit (can depend on x)
            #epsabs=1e-2,   # <-- Increased absolute tolerance
            #epsrel=1e-2,    # <-- Increased relative tolerance
        )[0] \
        + v[t_step] ** 2 / 2

    obj_value = np.sum(obj) * T / t_division - (obj[0] + obj[-1]) * T / t_division / 2
    print(obj_value)
    v = np.maximum(v - v_grad,0)
    print(v_grad)
    print(v)