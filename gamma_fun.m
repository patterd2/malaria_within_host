function y = gamma_fun(x)
    % delay function for bursting of red blood cells
    a1 = 8; % shape parameter
    b1 = 1/4; % scale parameter, mean should be alpha
    % mean = a1*b1
    y = gampdf(x,a1,b1);
end
