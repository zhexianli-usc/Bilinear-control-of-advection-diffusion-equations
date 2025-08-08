
T = 4;
L_1 = 2;
L = 10;

t_step = 20;
x_step = 20;
x_step_con = 2;


dt = T / t_step; 
dx = L_1 / x_step;
dx_con = L_1 / x_step_con;  
x = [1]; %[0, 0.5, 1];

v_opt = [zeros(1,6) zeros(1,4) zeros(1,5) zeros(1,5)];
state_v = zeros(t_step, x_step);
for i = 1 : t_step
    for j = 1 : x_step
        state_v(i, j) = integral (@(r) state_real(r,j * dx,dt * i,v_opt(1:t_step),dt,t_step),0,L);
    end
end
  

function obj = obj_deri(dt, var ,step, t_step, L)
    obj = var(step); %+ quad2d(@(x,r) derivative_real(r,x,dt * step,var(1:t_step),var(step),dt,t_step), 0, L, 0, L);

end

function bc = boundary_condition(t) 
    bc = 2 * abs(sin( pi * t / 2));
    % bc_1 = 5;
    % bc_2 = 0.5;
    % if t < 1/2
    %     bc = bc_1;
    % elseif t < 1
    %     bc = bc_2;
    % elseif t < 3/2
    %     bc = bc_1 - 1;
    % else
    %     bc = bc_2;
    % end
end

function matrix = necessary_stateub(var, x, dt, t_step,L)
    constraint_time = t_step / 2;
    matrix = zeros(3 * t_step,1);
    state_max = zeros(t_step,1);
    state_max(1:constraint_time) = 0.63; % sin(dt * (1:constraint_time)) * 0.4;
    state_max(constraint_time:t_step) = 0.63;
    v_min = 0;
    len_x = length(x);

    for i = 1 : t_step
        for j = 1 : len_x
            matrix(i) = var(i + j * t_step) * obj_deri(dt, var, i, t_step, L); 
            matrix(i + j * t_step) = (integral(@(r) state_real(r,x(j),dt * i,var(1:t_step),dt,t_step),0,L) + var(i + j * t_step)^2 - state_max(i));
        end
        matrix(i + (len_x + 1) * t_step) = (var(i) - var(i + (len_x + 1) * t_step)^2 - v_min);
    end
end

function matrix = necessary_statelb(var, x, dt, t_step,L)
    constraint_time = t_step / 2;
    matrix = zeros(3 * t_step,1);
    state_min = zeros(t_step,1);
    state_min(1:constraint_time) = sin(dt * (1:constraint_time)) * 0.4;
    state_min(constraint_time:t_step) = 0.35;
    v_max = 1;
    len_x = length(x);

    for i = 1 : t_step
        for j = 1 : len_x
            matrix(i) = var(i + j * t_step) * obj_deri(dt, var, i, t_step, L); 
            matrix(i + j * t_step) = 10 * (integral(@(r) state_real(r,x(j),dt * i,var(1:t_step),dt,t_step),0,L) - var(i + j * t_step)^2 - state_min(i));
        end
        matrix(i + (len_x + 1) * t_step) = 10 * (var(i) + var(i + (len_x + 1) * t_step)^2 - v_max);
    end
end
% quad2d(@(x,r) realpart(r,x,t,v,v(1),dt,t_step), 0, L, 0, L)
% integral(@(r) realpart(r,x,t,v,v(1),dt,t_step), 0, L)

function integrand = derivative_real(k,x,t,v,v_cur,dt,t_step)
    r1 = exp(1i*pi/8);
    r2 = exp(7i*pi/8);

    integrand = real(Vv(k * r1,x,t,v,v_cur,dt,t_step) * r1 - Vv(k * r2,x,t,v,v_cur,dt,t_step) * r2) / (2*pi);
                % .* 
                % (- real(V(k * r1,x,t,v,dt,t_step) * r1 - V(k * r2,x,t,v,dt,t_step) * r2) /(2*pi));
end

function integrand = state_real(k,x,t,v,dt,t_step)
    r1 = exp(1i*pi/8);
    r2 = exp(7i*pi/8);
    integrand = - real(V(k * r1,x,t,v,dt,t_step) * r1 - V(k * r2,x,t,v,dt,t_step) * r2) /(2*pi);
end

function state = V(k,x,t,v,dt,t_step)
    v_tilde = vint(t,v,dt,t_step);
    state   = exp(1i * k .* x ) .* (2i * k - v_tilde / t) .* B(k,v,t,dt,t_step);
end

function derivative = Vv(k,x,t,v,v_cur,dt,t_step)
    v_tilde    = vint(t,v,dt,t_step);
    derivative = exp(1i * k .* x ) .* ((2 * k.^2 * v_cur + 1i * k * v_cur * v_tilde / t) + (v_cur / t))...
        .* B(k,v,t,dt,t_step) ...
        + exp(1i * k .* x ) .* (2 * k.^2 + 1i * k * v_tilde/t ) .* Bv(k,v,t,dt,t_step);
end

function g0 = B(k,v,t,dt,t_step)
    dtt = 0.01;
    step = floor(t/dtt);
    g0 = 0;
    for i = 1 : step
        tt = dtt * i;
        if i == 1 || step
            bound = 1/2;
        else
            bound = 1;
        end
        g0 = g0 + bound * boundary_condition(tt) * dtt * exp(k.^2 * (tt-t) + 1i * vint(tt,v,dt,t_step) * k);
    end
    g0 = g0 .* exp(- 1i * vint(t,v,dt,t_step) * k);
    g0 = g0 + (t - tt) * boundary_condition(t);
end

function g0 = Bv(k,v,t,dt,t_step)
    dtt = 0.01;
    step = floor(t/dtt);
    g0 = 0;
    for i = 1 : step
        tt = dtt * i;
        if i == 1 || step
            bound = 1/2;
        else
            bound = 1;
        end

        v_step = min(floor(tt / dt) + 1, t_step);
        g0 = g0 + bound * boundary_condition(tt) * v(v_step) * dtt * exp(k.^2 * (tt-t) + 1i * vint(tt,v,dt,t_step) * k);
    end
    g0 = g0 .* exp(- 1i * vint(t,v,dt,t_step) * k);
    v_step = min(floor(tt / dt) + 1, t_step);
    g0 = g0 + v(v_step) * (t - tt) * boundary_condition(t);
end


function vi = vint(t,v,dt,t_step)
    step = floor(t / dt);
    vi = sum(v(1:step) * dt);
    if step < t_step
        vi = vi + (t - step * dt) * v(step + 1);
    end
end