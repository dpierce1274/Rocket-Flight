function [rho] = Density(h)

global rhos;

rho = rhos*exp((-1/8000)*h); % Function for changing atmospheric density

end

