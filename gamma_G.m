function y = gamma_G(x)
    % delay function for bursting of red blood cells
    a1 = 5; % shape parameter
    b1 = 1/4.47; % scale parameter
    % mean = a1*b1
    y = gampdf(x,a1,b1); 
end
