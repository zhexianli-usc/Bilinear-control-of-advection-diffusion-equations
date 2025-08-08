# Bilinear control of advection diffusion equations

This repository includes codes for computing constrained bilinear control of advection diffusion equations.
The folder "Matlab" includes all MATLAB files. The main file for the code is "IntegralEquation.m"
The folder "Figure-data" includes data used in plots.

Consider the following 1-d advection-diffusion-reaction equation on the domain $0\leq t\leq T,x\geq0$:\
$$\frac{\partial}{\partial t}\psi(x,t) + v \frac{\partial}{\partial x}\psi(x,t) = D \frac{\partial^2}{\partial x^2}\psi(x,t) - u(t) \psi(x,t)$$\
subject to initial and boundary conditions $\psi(x,0)=\phi_o(x),\psi(0,t)=h(t),\psi(\infty,t)=0$. 

Our goal is to design $u(t)$ such that certain objective is minimized.
An example of objective is the following:\
$$J(u)=\int_{0}^{T}w_1\left[\int_{0}^{\infty}\psi(x,t) dx\right]  + w_2u^2(t)dt$$
