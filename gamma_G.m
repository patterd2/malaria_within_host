function y = gamma_G(x)
    % delay function for bursting of red blood cells
    a1 = 24*33; % shape parameter
    b1 = 1/3; % scale parameter, mean should be alpha_G
    % mean = a1*b1
    %x = linspace(0,20*24,1000);
    y = gampdf(x,a1,b1); 
    %plot(x/24,y)
end
