function y = phi(x,IT,s)
% delay function for bursting of red blood cells
%phi0 = 1./(1+exp(-(0-IT)./s));
%y = 1./(1+exp(-(x-IT)./s)) - phi0; % functional form taken from REF? need to update
y = 0.5*(1 + tanh((x-IT)/s));
end
