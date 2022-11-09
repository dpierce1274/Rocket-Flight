function [ T ] = Thrust( t )
% The Thrust Function determines the thrust of the rocket at any given time
             
global array tend tdata
tdata = array(1:end,1);   % Time Data from excel spreadsheet (s)
Tdata = array(1:end,2);   % Thrust Data from excel spreadsheet (N)
tend = array(end,1);

if t < tend
    T = interp1(tdata,Tdata,t,'linear');    % Interpolation to get thrust (N) at each time increment        
else 
    T = 0;
end

