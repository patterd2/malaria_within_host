function y = gamma_fun(x,h)
    % delay function for bursting of red blood cells
    a1 = 24*8; % shape parameter
    b1 = 1/4; % scale parameter, mean should be alpha
    % mean = a1*b1
    pi1 = 1-gamcdf(x,a1,b1);%1000*gampdf(x,a1,b1);
    pi2 = 1-gamcdf(x+h,a1,b1);%1000*gampdf(x+h,a1,b1);

    y = 1/h*log(pi1./pi2);
    id = find(y>1,1,'first');
    y(id:end) = 1;
end
