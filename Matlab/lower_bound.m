

% v = flip([1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0]);
lam = [0    0.05    0.1 0.15 0.2 0.25 0.3];
ub = [0.7709 0.5240 0.363 0.2562 0.1844 0.1352 0.1009];
lam_solute = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
ub_solute = [1.5541 1.1203 0.8714 0.7039 0.5867 0.4981 0.4277 0.3705 0.3235 0.2844 0.2514];
% bound = 2 * flip([0.382 0.371 0.36 0.348 0.336 0.324 0.311 0.299 0.286 0.273 0.260]);
% 
% upper = flip([0.6789 0.6803 0.6850 0.6927 0.7066 0.7139 0.7192 0.7256 0.7275 0.7347 0.7417]);
% 
x = [0 dt * (1:t_step)];
bc = zeros(1, t_step + 1);
constraint_time = t_step / 2;
state_min = zeros(1,t_step + 1);
state_min(1:constraint_time+1) = 0.4 * exp(0.5 - [0 dt * (1:constraint_time)]);
state_min(constraint_time+1:t_step+1) = 0.4;
hfig = figure;

% 2) Create your figure and plot


% plot(x, state_min, 'LineWidth', 3, 'Marker', 'o');
% plot(lam, ub, 'LineWidth', 4,'Color',[0.9290, 0.6940, 0.1250],'Marker','o');
% % plot(v, upper, 'LineWidth', 3, 'Color',[0.9290, 0.6940, 0.1250],'Marker','o');

%%% Plotting boundary condition curve

% for i = 1 : t_step + 1
%     bc(i) = boundary_condition(x(i));
% end
% plot(x, bc, 'LineWidth', 4, 'Marker','o');


%%% Plotting comparison of solutions when fixing x = 1  %%%%%%

% plot(dt * (0:t_step), [state(1:end)'] , 'LineWidth',4)
% hold on
% plot(0.01 * (0:50:1000), sol(1:50:1001, 101)' , 'Marker','o','MarkerSize',12, 'LineStyle', 'none','LineWidth',3)
% hold on
% legend( 'Unified Transform','Benchmark','Location','best','Interpreter','latex', ...
%         'FontSize', 21)

%%% Plotting comparison of solutions when fixing t = 1  %%%%%%
plot(dx * (1:x_step), [state(1:end)'] , 'LineWidth',4)
hold on
plot(0.01 * (10:50:1000), sol(101, 10:50:1000),'Marker','o','MarkerSize',12, 'LineStyle', 'none','LineWidth',3)
hold on
legend( 'Unified Transform','Benchmark','Location','best','Interpreter','latex', ...
        'FontSize', 21)


% % 3) Customize labels and title with font name & size
xlabel('$x$', ...
        'Interpreter','latex', ...
        'FontSize', 25);
% 'Upper bound on velocity $\bar{\textit{v}}$'
% 'Boundary condition $\phi(0,t)$'
% '$\phi(1,t)$'
        ylabel('$C(x,1)$' ,...
                    'Interpreter','latex', ...
        'FontSize', 25);
            % 'Maximally feasible $\mathit{\mathbf{\phi}}_{\max}$', ...
% title('Damped Sinusoid and Cosinusoid', ...
%       'FontSize', 18, ...
%       'FontName', 'Helvetica', ...
%       'FontWeight', 'bold');
% 
% 4) Adjust the axes tick labels
ax = gca;
ax.FontName   = 'Arial';         % tick-label font
ax.FontSize   = 25;
ax.Box        = 'on';            % draw box around axes
ax.TickLabelInterpreter = 'latex';

% 5) Optional: set global defaults for all future plots in this session
set(groot, ...
    'defaultAxesFontName',   'Arial')
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 25)
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'Latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(findall(hfig, '-property', 'LindWidth'), 'LindWidth', 4)

% Now any new figure you open will inherit these settings!

exportgraphics(gca, '/Users/zhexianli/Dropbox/Apps/Overleaf/Control in mechanics/fig/solution-compare-solute-x-1.pdf', ...
    'BackgroundColor','none', ...
    'ContentType','vector');

% % % 
% exportgraphics(gca, '/Users/zhexianli/Dropbox/Apps/Overleaf/Control in mechanics/fig/maximal-upper-bound.pdf', ...
%     'BackgroundColor','none', ...
%     'ContentType','vector');
% % 
% exportgraphics(gca, '/Users/zhexianli/Dropbox/Apps/Overleaf/Control in mechanics/fig/Boundary-condition/boundary-condition-nuclear-sin.pdf', ...
%     'BackgroundColor','none', ...
%     'ContentType','vector');