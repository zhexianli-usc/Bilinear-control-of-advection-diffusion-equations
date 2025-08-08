function bc = boundary_condition(t) 
    % len_t = length(t);
    % bc = -ones(1,len_t);
    % bc = besselj(0, t) + 1;
    % bc = 2 * abs(sin(2 * pi * t));
    bc = - 10 * (1 + sin(2 * pi * t) / 3);
    % bc = - 10;
    % bc = 1 + cos( pi * (t)) / 3;
    % if t < 1.5
    %     bc = 2 * abs(sin(pi * t /2));
    % elseif t < 2
    %     bc = abs(sin(pi * (t+1) /4)) ;
    % elseif t < 2.5
    %     bc = abs(sin(pi * (t-1) /4)) ;
    % else
    %     bc = 2 * abs(sin(pi * t /2));
    % end
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