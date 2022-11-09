function [ dvdt ] = acceleration1(tb,v,h)
% The acceleration function is used to calculate the acceleration of a
% rocket at a given time and velocity.

global g0;
 
dvdt = (Thrust(tb) - Drag(v,h))/Mass1(tb) - g0;



end