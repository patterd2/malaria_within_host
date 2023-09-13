function y = phi(x,IT,s)
% delay function for bursting of red blood cells
   y = 25./(1+exp(x-IT)./s); % functional form taken from REF? need to update
end
