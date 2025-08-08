% Mesh and time span
xmesh = linspace(0,100,10001); % 100 points between x=0 and x=1
tspan = linspace(0,100,10001); % Solve for 50 time steps up to t=1

% Solve the PDE
sol = pdepe(0, @pdefun, @icfun, @bcfun, xmesh, tspan);

% Visualize the solution
% figure;
% surf(xmesh, tspan, sol)
% [C, h] = contourf(xmesh, tspan, sol);
% xlabel('x');
% ylabel('t');
% zlabel('u');
% title('Heat Equation Solution');

% Define the PDE function (pdefun)
function [c,f,s] = pdefun(x,t,u,DuDx)
  alpha = 1; % Thermal diffusivity
  c = 1;
  f = DuDx;
  % s = - alpha * 9.21 * 0.1532 * exp(-t) * u; % No source term
  s = - DuDx - 0.001 * 11.4 / 0.426^2 * (1 + sin(t) / 2) * u;
end

% Define the initial condition (icfun)
function u0 = icfun(x)
  u0 = 2 * exp(-1 * x); % Example initial temperature
end

% Define the boundary conditions (bcfun)
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
  pl = ul - (besselj(0, t) + 1);
  ql = 0;
  % pl = 10 * (1 + sin(t)/2);
  % ql = 9.21; %(besselj(0, t) + 1);
  pr = ur;
  qr = 0; % Zero heat flux at the boundaries
end