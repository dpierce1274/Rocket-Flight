%% Modified Barrowman Stability Equation
% Solves Barrowman Stability equation for NPS High-Power Rocket

%% Program Start

clear, clc, format compact

%% Rocket Inputs

calc_twostage = false;   % Calculate sustainer alone or two stage configuration

% Sustainer Airframe and Fin Inputs

Ln = 34.5;              % Length of nosecone (in)
d = 6.17;               % Diameter at base of nose (in)
df = 6.17;              % Diameter at front of transition (in)
dr = 4.024;             % Diameter at rear of transition (in)
Lt = 4;                 % Length of transition (in)
Xp = 34.5;              % Distance from tip of nose to front of transition (in)
Cr = 10;                % Fin root chord (in)
Ct = 2;                 % Fin tip chord (in)
S = 4.625;              % Fin semispan (in)
Lf = 5.58;              % Length of fin mid-chord line (in)
R = 2.01;               % Radius of body at aft end (in)
Xr = 7.125;             % Distance between fin root leading edge and fin tip leading edge parallel to body (in)
Xb = 109.5;             % Distance from nose tip to fin root chord leading edge (in)
N = 3;                  % Number of fins

% Booster Airframe and Fin Inputs

b_d = 6.17;               % Diameter at base of nose (in)
b_df = 4.024;             % Diameter at front of transition (in)
b_dr = 6;                 % Diameter at rear of transition (in)
b_Lt = 3;                 % Length of transition (in)
b_Xp = 119.5;             % Distance from tip of nose to front of transition (in)
b_Cr = 13.5;              % Fin root chord (in)
b_Ct = 2.6;               % Fin tip chord (in)
b_S = 6.5;                % Fin semispan (in)
b_Lf = 7.73;              % Length of fin mid-chord line (in)
b_R = 3;                  % Radius of body at aft end (in)
b_Xr = 9.625;             % Distance between fin root leading edge and fin tip leading edge parallel to body (in)
b_Xb = 245.6;             % Distance from nose tip to fin root chord leading edge (in)
b_N = 3;                  % Number of fins

% Nosecone Parameter Inputs

Rad = 6.17/2;           % Radius of nosecone at base (in)
L = 34.5;               % Length of nosecone (in)
fidelity = 1000;        % Fidelity of calculation

%% Von Karman LD-Haack Nosecone Term Calculations

theta = zeros(1,fidelity);
y = zeros(1,fidelity);
x = linspace(0,L,fidelity);
neg_y = zeros(1,fidelity);
vol_y = zeros(1,fidelity);
dx = L/fidelity;
y(1) = 0;

for k = 1:fidelity
    theta(k) = acos(1-(2*x(k))/L);
    y(k) = (R/sqrt(pi))*sqrt(theta(k)-(sin(2*theta(k)))/2);
    neg_y(k) = -(R/sqrt(pi))*sqrt(theta(k)-(sin(2*theta(k)))/2);
    vol_y(k+1) = vol_y(k)+pi*y(k)^2*dx;

end

volume = vol_y(end);

Z = (L-(volume/(pi*y(end)^2)));

figure(1)
plot(x(k)*(Z/L),0,'or',x,y,x,neg_y)
title('LD-Haack (Von Karman) Nosecone'), xlabel('Length (in)'), ylabel('Radius (in)')
axis('equal')
legend('CP Location')
%% Nose Cone Term

Xn = (Z/L)*Ln;          % For 
Cnn = 2;

%% Conical Transition Terms
% Sustainer Terms
Cnt = 2*((dr/d)^2-(df/d)^2);
Xt = Xp+(Lt/3)*(1+((1-(df/dr))/(1-(df/dr)^2)));

% Booster Terms
b_Cnt = 2*((b_dr/b_d)^2-(b_df/b_d)^2);
b_Xt = b_Xp+(b_Lt/3)*(1+((1-(b_df/b_dr))/(1-(b_df/b_dr)^2)));

%% Fin Terms

Cnf = (1+(R/(S+R)))*((4*N*(S/d)^2)/(1+(sqrt(1+(((2*Lf)/(Cr+Ct))^2)))));
Xf = Xb+(((Xr/3)*(Cr+2*Ct))/(Cr+Ct))+((1/6)*((Cr+Ct)-((Cr*Ct)/(Cr+Ct))));

b_Cnf = (1+(b_R/(b_S+b_R)))*((4*b_N*(b_S/b_d)^2)/(1+(sqrt(1+(((2*b_Lf)/(b_Cr+b_Ct))^2)))));
b_Xf = b_Xb+(((b_Xr/3)*(b_Cr+2*b_Ct))/(b_Cr+b_Ct))+((1/6)*((b_Cr+b_Ct)-((b_Cr*b_Ct)/(b_Cr+b_Ct))));

%% Center of Pressure

if calc_twostage == true
    Cnr = Cnn + Cnf + Cnt + b_Cnf + b_Cnt;
    X = ((Cnn*Xn)+(Cnf*Xf)+(Cnt*Xt)+(b_Cnf*b_Xf)+(b_Cnt*b_Xt))/Cnr; 
else
    Cnr = Cnn + Cnf + Cnt;
    X = ((Cnn*Xn)+(Cnf*Xf)+(Cnt*Xt))/Cnr;               % Location from Nose Tip
end


fprintf('LD-Haack Barrowman Parameter = %.3f \n',Z/L)
fprintf('Center of Pressure Distance from Nose Cone = %3.2f in\n', X)