function [ m ] = Mass( t )
%The Mass function gives the mass of the rocket at a given time after
%launch. A linear relationship with time.

global m1 motor pro tb tend

mint = m1+motor+pro;         % Initial Mass (Empty Rocket + Empty Engine + Propellant)
mfin = m1+motor;             % Final Mass (Empty Rocket + Empty Engine)
dmdt = (mint-mfin)/tb(end);  % Linear Mass change over burnout time

if t <= tend                 % Losing Mass during burn phase
    m = -t.*dmdt+mint;
else
    m = mfin;                % Consant Mass during coast phase
end

end

