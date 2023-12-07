function y = gammaV(x,h)
    % infectivity function
    a1 = 24*47; % shape parameter
    b1 = 1/4.47; % scale parameter
    % mean = a1*b1
    pi1 = 1-gamcdf(x,a1,b1);
    pi2 = 1-gamcdf(x+h,a1,b1);

    y = (1/h)*log(pi1./pi2);
    id = find(y>.8,1,'first');
    y(id:end) = .8;
end
