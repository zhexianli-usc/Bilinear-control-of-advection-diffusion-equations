function ic_hat = initial_condition(k) 
    ic_hat =  1 ./ (1i * k + 10); % for phi_0(x) = exp(-x)
    % ic_hat = 0; % for phi_0(x) = 0
    % ic_hat = 1; % for phi_0(x) = 1
end