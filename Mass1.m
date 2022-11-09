function [ m ] = Mass1( t )
%The Mass function gives the mass of the rocket at a given time after
%launch. A linear relationship with time.

global m1 O8000 pro tb tend

mint1 = m1+O8000*2+pro*2;         % Initial Mass (Empty Rocket + Empty Engine + Propellant)
mfin1 = m1+O8000*2+pro;
mint2 = m1+O8000+pro;             % Final Mass (Empty Rocket + Empty Engine)
mfin2 = m1+O8000;
dmdt1 = (mint1-mfin1)/tb(end);    % Linear Mass change over burnout time
dmdt2 = (mint1-mfin1)/tb(end);

if t <= tend                 % Losing Mass during burn phase
    m = -t.*dmdt1+mint1;
elseif t >= tend + 0.5
    m = mint2;
elseif t >= tend + 4
    m = -t.*dmdt2+mint2;
elseif t >= tend*2 + 4
    m = mfin2;                % Consant Mass during coast phase
end

end