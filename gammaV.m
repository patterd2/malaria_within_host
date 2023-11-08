function y = gammaV(x)
    % infectivity function
    a1 = 24*47; % shape parameter
    b1 = 1/4.47; % scale parameter
    % mean = a1*b1
    y = gampdf(x,a1,b1); 
end
