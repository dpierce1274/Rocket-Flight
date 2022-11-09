function [rho] = Rho(h)

global rhos

rho(k) = rhos*exp((-1/8000)*h);

end

