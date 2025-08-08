% 1) Define your grid
t = dt * (1:t_step);
x = dx * (1:x_step);
[T, X] = meshgrid(t, x);

% 2) Define your function Z = f(X,Y)

% 3) Create the contour plot
hfig = figure;
% contourf for filled contours, contour for lines only
[C, h] = contourf(x, t, state, ...      % data
                  1000, ...           % number of contour levels
                  'LineColor','none'); % no black lines between levels

% 4) Add a colorbar
cb = colorbar;
cb.Label.String = '$C(x,t)$';
cb.Label.Interpreter = 'latex';    % if you want LaTeX

% 5) Customize the colormap (optional)
colormap(hot);  % try parula, jet, viridis (if you have it), etc.

% 6) Label axes and title
xlabel('$x$', ...
       'Interpreter','latex', ...
       'FontSize',25);
ylabel('$t$', ...
       'Interpreter','latex', ...
       'FontSize',25);
       % title('Filled Contour of $\phi(t,x)$', ...
%       'Interpreter','latex', ...
%       'FontSize',18);

% 7) Adjust axes properties for better readability
ax = gca;
ax.TickLabelInterpreter = 'latex';


% 8) (Optional) Add contour labels on lines
% If you prefer line contours and numeric labels:
hold on;
[~, hLine] = contour(x, t, state, [.5 .5], 'LineColor','w', 'LineWidth', 3); 
clabel(C, hLine, 'FontSize',25, 'LabelSpacing',300, 'Color','w', 'Interpreter', 'Latex');

hold on;

% plot(ones(1,t_step),t(1:t_step), 'y--', 'LineWidth',3);
% plot(1.99 * ones(1,t_step),t(1:t_step), 'y--', 'LineWidth',3);
% plot(x(1:t_step),3.99 * ones(1,t_step), 'y--', 'LineWidth',3);


% hold on;
% plot(x(1:x_step/2),ones(1,x_step/2), 'w--', 'LineWidth',2);
% hold on;
% plot(x(1:x_step/2),ones(1,x_step/2) + 0.6, 'w--', 'LineWidth',2);


% plot(x(x_step/2:x_step),ones(1,x_step/2+1) / 2 , 'w--', 'LineWidth',2)


set(findall(hfig, '-property', 'FontSize'), 'FontSize', 25);
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'Latex');
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
set(findall(hfig, '-property', 'LindWidth'), 'LindWidth', 3);
hold off;
% 
% exportgraphics(gca, '/Users/zhexianli/Dropbox/Apps/Overleaf/Control in mechanics/fig/contour/state-contour-nuclear-constant.pdf', ...
%     'BackgroundColor','none', ...
%     'ContentType','vector')

