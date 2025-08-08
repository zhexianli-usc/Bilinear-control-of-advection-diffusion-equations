
T = 1;
L_1 = 1;
L = 10;

t_step = 10;
x_step = 100;
x_step_con = 100;


dt = T / t_step; 
dx = L_1 / x_step;
dx_con = L_1 / x_step_con;  
x = [0.5]; %[0, 0.5, 1]; /

lam0  = zeros((length(x) + 2) * t_step,1) + 1;
options=optimset('Display','off','OutputFcn', @printX_perIter,'MaxIter',100000,'MaxFunEvals',100000,'TolFun',0.0001, 'Algorithm','trust-region-reflective');
% fun = @(v) necessary_statelb(v,x, dt, t_step,L);
fun = @(lam) necessary_stateub(lam,x, dt, t_step,L,T);


[lam,fval] = fsolve(fun, lam0, options);
 
lam_opt = lam(1:t_step);
% lam_opt = ones(1,t_step) * 0.0 * D;
% v_opt = [linspace(0.1, 2, 20) ones(1, t_step - 20)];
% v_opt = abs(sin(linspace(dt, T, t_step))) / 2;
% lam_opt = 0.1523 * 9.21 * (1 + sin(linspace(0, T, t_step + 1)) / 2);
% v_opt   = zeros(1, t_step);
% state_x = 1;

% state = zeros(t_step, x_step);
state = zeros(t_step + 1, 1);

for i = 1 : t_step + 1
    % for j = 1 : x_step
        state(i) = integral_rep(x, dt * (i - 1), lam_opt(1:t_step), dt, t_step, 10 * L);
        % state(i)
        % state(j) = integral_rep(dx * (j), 5, lam_opt(1:t_step), dt, t_step, L);
        % state(i, j) = integral_rep(dx * j, dt * i, lam_opt(1:t_step), dt, t_step, L);
    % end
end

% state = zeros(x_step, 1);
% state_t = 1;
% for j = 1 : x_step
%     state(j) = integral (@(r) state_real(r, dx * j, state_t, v_opt(1:t_step), dt, t_step),0,100);
% end
  
function stop = printX_perIter(x, optimValues, state)
    stop = false;  % returning true would halt fsolve
    if strcmp(state,'iter')
        % Print iteration number + full x-vector
        fprintf('Iter %2d: x = [%s], ||F|| = %.3e\n', ...
                optimValues.iteration, ...
                num2str(x','%.6g '), ...
                norm(optimValues.fval));
    end
end

function obj = obj_deri(x, dt, lam, Lagrange_state, Lagrange_con, step, t_step, L, T)
    G_u = 0;
    tau_m = dt * step;
    for i = step : t_step
        tau_i = dt * i;
        G_u = G_u + integral2(@(r, t) - Lagrange_state(i) * der_integral(r, x, t, lam, dt, t_step, tau_m), 0, 10 * L, tau_i, tau_i + dt);
    end
    obj = lam(step) + G_u;% - Lagrange_con(step);

end

function derivative = Vv(k,x,t,lam,dt,t_step,tau)
    v_tilde = vint(t);
    D = 1 / 9.21 ;
    lam_tilde = lam_int(t,lam,dt,t_step);
    omega = D * k.^2 .* t + 1i * v_tilde * k + lam_tilde;

    derivative   = exp(1i * k * x - omega) .* initial_condition(k) ...
              - exp(1i * k .* x ) .* ( exp(-omega) .* -initial_condition(-k) ...
                                     + ( 2 * D .* B(k, lam, D, t, tau, dt, t_step, lam_tilde)));
    % lam_tilde = lam_int(t,lam,dt,t_step);
    % derivative = exp(1i * k * x - k.^2 .* t - 1i * v_tilde .* k - lam_tilde) .* initial_condition(k) ...
    %           - exp(1i * k * x ) .* ( exp(-k.^2 .* t - 1i * v_tilde .* k - lam_tilde) .* initial_condition(-k - 1i) ...
    %                                  + (2i * k - 1) .* B(k, lam, t, tau, dt, t_step, lam_tilde));
end

function integrand = der_integral(k,x,t,lam,dt,t_step,tau) 
    r1 = exp(1i*pi/8);
    r2 = exp(7i*pi/8);
    integrand = real(Vv(k * r1,x,t,lam,dt,t_step,tau) * r1 - Vv(k * r2,x,t,lam,dt,t_step,tau) * r2) / (2 * pi);
end

function matrix = necessary_stateub(var, x, dt, t_step, L, T)
    len_x = length(x);
    constraint_time = t_step / 2;
    matrix = zeros((len_x + 2) * t_step,1);
    state_max = zeros(t_step,1);
    state_max(1:constraint_time) = 0.5;% * exp(x(end) - dt * (1:constraint_time)); % sin(dt * (1:constraint_time)) * 0.4;
    state_max(constraint_time:t_step) = 0.5;
    v_min = 0;
    v_max = 3;
    % var contains control u = var(1:t_step), 
    % multiplier for state constraint lagrange = var(t_step + 1 : 2 * t_step)
    % multiplier for control constraint lagrange = var(2 * t_step + 1 : 3 * t_step)
    % auxilary variable for state constraint = var(3 * t_step + 1 : 4 * t_step)
    % auxilary variable for control constraint = var(4 * t_step + 1 : 5 * t_step)
    
    matrix = zeros((len_x + 2) * t_step,1);
    lam            = var(1 : t_step);
    Lagrange_state = var(t_step + 1 : 2 * t_step);
    aux_state      = var(2 * t_step + 1 : 3 * t_step);
    % Lagrange_con   = var(3 * t_step + 1 : 4 * t_step);
    % aux_con = var(4 * t_step + 1 : 5 * t_step);
    for i = 1 : t_step
        for j = 1 : len_x
            % Derivative condition
            matrix(i) = obj_deri(x(j), dt, lam, Lagrange_state, [], i, t_step, L, T) ;
            % state constraint
            matrix(i + j * t_step) = integral_rep(x(j),dt * i, lam, dt,t_step,L) + aux_state(i)^2 - state_max(i);
        end
        matrix(i + (len_x + 1) * t_step) = aux_state(i) * Lagrange_state(i);
        %   control constraint
        % matrix(i + (len_x + 2) * t_step) = lam(i) - aux_con(i)^2 - v_min;
        % matrix(i + (len_x + 3) * t_step) = aux_con(i) * Lagrange_con(i);
        
        % matrix(i + (len_x + 2) * t_step) = (var(i) + var(i + (len_x + 2) * t_step)^2 - v_max);
    end
end

function matrix = necessary_statelb(var, x, dt, t_step,L)
    constraint_time = t_step / 2;
    matrix = zeros(3 * t_step,1);
    state_min = zeros(t_step,1);
    % state_min(1:constraint_time) = sin(dt * (1:constraint_time)) * 0.7;
    state_min(constraint_time:constraint_time + 1) = 0.6;
    v_max = 1;
    v_min = 0;
    len_x = length(x);

    for i = 1 : t_step
        for j = 1 : len_x
            matrix(i) = var(i + j * t_step) * obj_deri(dt, var, i, t_step, L); 
            if i == constraint_time || i == constraint_time + 1
                matrix(i + j * t_step) =  integral_rep(x(j), dt * i, var(1:t_step), dt, t_step, L) - var(i + j * t_step)^2 - state_min(i);
            end
        end
        matrix(i + (len_x + 1) * t_step) =  (var(i) + var(i + (len_x + 1) * t_step)^2 - v_max);
          % matrix(i + (len_x + 1) * t_step) =  (var(i) - var(i + (len_x + 1) * t_step)^2 - v_min);
    end
end
% quad2d(@(x,r) realpart(r,x,t,v,v(1),dt,t_step), 0, L, 0, L)
% integral(@(r) realpart(r,x,t,v,v(1),dt,t_step), 0, L)

function state_val = integral_rep(x, t, lam, dt, t_step, L)
    state_val = integral(@(r) state_real(r,x, t,lam,dt, t_step),0, L);%...
                % + integral(@(k) exp(1i * k * x - k.^2 * t - 1i * vint(t,v,dt,t_step) * k) .* initial_condition(k), -10 * L, 10 * L);
end

function integrand = derivative_real(k,x,t,v,v_cur,dt,t_step)
    r1 = exp(1i*pi/8);
    r2 = exp(7i*pi/8);

    integrand = real(Vv(k * r1,x,t,v,v_cur,dt,t_step) * r1 - Vv(k * r2,x,t,v,v_cur,dt,t_step) * r2) / (2*pi);
                % .* 
                % (- real(V(k * r1,x,t,v,dt,t_step) * r1 - V(k * r2,x,t,v,dt,t_step) * r2) /(2*pi));
end

function integrand = state_real(k,x,t,lam,dt,t_step)
    r1 = exp(1i*pi/8);
    r2 = exp(7i*pi/8);
    integrand = real(V(k * r1,x,t,lam, dt,t_step) * r1 - V(k * r2,x,t,lam, dt,t_step) * r2) /(2*pi);
end

function state = V(k,x,t,lam,dt,t_step)
    v_tilde = vint(t);
    D = 1 / 9.21 ;
    lam_tilde = lam_int(t,lam,dt,t_step);
    omega = D * k.^2 * t + 1i * v_tilde * k + lam_tilde;

    state   = exp(1i * k * x - omega) .* initial_condition(k) ...
              - exp(1i * k .* x ) .* ( exp(-omega) .* -initial_condition(-k) ...
                                     + ( 2 * D .* B(k,lam,D,t,t,dt,t_step, lam_tilde)));
end

function g0 = B(k, lam, D, t, tau, dt, t_step, lam_tilde)
    % g0 = (1- exp(-k.^2 * t)) ./ k.^2;
    % len_k =  length(k);
    % g0 = zeros(1, len_k);
    % for i = 1 : len_k
    %     g0(i) = integral(@(tt) boundary_condition(tt) ...
    %                             .* exp(lam_int(tt,lam,dt,t_step) + D * k(i)^2 * (tt - t) ...
    %                              + 1i * (vint(tt) - vint(t)) * k(i)),0,t,'RelTol',1e-9);
    % end
    % g0 = g0 * exp(-lam_int(t,lam,dt,t_step));
    g0 = zeros(size(t));
    if tau > 0 
        dtt = 0.01;
        step = floor(tau / dtt);
        for i = 1 : step
            tt = dtt * i;
            if i == 1 || i == step
                bound = 1 / 2;
            else
                bound = 1;
            end
            g0 = g0 + bound * boundary_condition(tt) * dtt ...
                            * exp(lam_int(tt,lam,dt,t_step) - lam_tilde ...
                            + D * k.^2 .* (tt-t) + 1i * (vint(tt) - vint(t)) .* k);
        end
        g0 = g0 + (tau - tt) * (boundary_condition(tt) ...
                            * exp(lam_int(tt,lam,dt,t_step) - lam_tilde + D * k.^2 .* (tt-t) ...
                                  + 1i * (vint(tt) - vint(t)) .* k)...
                               + boundary_condition(tau)) / 2;
    end
end

function lam_total = lam_int(t,lam,dt,t_step)
    % lam_total = 0;
    % lam_total = 9.21 * 0.1532 * (1-exp(-t));
    lam_total = zeros(size(t));

    for i = 1 : size(t,1)
        for j = 1 : size(t,2)
            if t(i, j) > 0
                step = floor(t(i, j) / dt);
                lam_total(i, j) = sum(lam(1:step) * dt);
                if step < t_step
                    lam_total(i, j) = lam_total(i, j) + (t(i, j) - step * dt) * lam(step + 1);
                end
            end
        end
    end
end

function vi = vint(t)
    vi = 0;
end