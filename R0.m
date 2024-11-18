baseline_parameter_set;

X_max = 1000*24; % max time in days
h = 0.001; % time/age step size in hours, same across all timescales

x = (0:h:X_max)';
nx = length(x);

gamma_int = h*cumsum(gamma_fun(x,h));

%P.sigma = 0.5/24;
fun = ((P.beta*(1-P.c).*P.p*P.Bbar)/(P.muM + P.p*P.Bbar)).*...
    gamma_fun(x,h).*exp(-(P.mu+P.sigma).*x - gamma_int );

%plot(x,fun);

disp(['R_0 = ',num2str(h*sum(fun))]);