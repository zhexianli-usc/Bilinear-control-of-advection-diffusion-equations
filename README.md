# Bilinear control of advection diffusion equations

This repository includes codes for implementing bilinear control of advection diffusion equations.

Consider the following 1-d advection diffusion equation on the domain $0\leq t\leq T,0<x<\infty$:
$$
\frac{\partial}{\partial t}\phi(x,t) + v(t) \frac{\partial}{\partial x}\phi(x,t) = \frac{\partial^2}{\partial x^2}\phi(x,t)
$$
subject to initial and boundary conditions $\phi(x,0)=0,\phi(0,t)=1,\phi(\infty,t)=0$. This situation represents, for example, advective-dispersive transport of a solute through a homogeneous, fully saturated, porous medium where an inactive tracer is continuously
injected at the origin with inlet concentration $1$, i.e., these conditions are met in laboratory scale tracer experiments in soil columns

Our goal is to design the function $v(t)$ such that certain objective is optimized. An example of objective is the following:
$$
J(v)=\int_{0}^{T}\left[\int_{0}^{\infty}\phi(x,t) dx\right]  - v^2(t)dt
$$
Maximizing $J(v)$ over $v$ can be interpreted as finding the minimum energy of $v$ to accelerate solute transport. 

The python file implemented a gradient descent method to solve optimal $v$. The gradient is computed in closed form using the Fokas method.