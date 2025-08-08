plot(0.1 * (0:1000), sol(1:end, 1000),  'LineWidth',3)

t_all = 0.1 * (0:1000);

x = -90;

heat_state = integral(1 / (sqrt(4 * pi * t_all)));