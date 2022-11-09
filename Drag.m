function [ D ] = Drag( v,h )
%The Drag function calculates the drag on an object at a certain velocity

global Cd A         % Calls the values for density, coefficient of drag, and area from the main program

D = .5*Density(h)*v^2*Cd*A;   % Function for drag (N)


end



