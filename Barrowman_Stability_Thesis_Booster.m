%% Computing Center of Pressure (Level 3)

clear, clc, format compact

%% Inputs

Ln = 2;                 % Length of nosecone (in)
d = 6.00;               % Diameter at base of nose (in)
df = 4.00;              % Diameter at front of transition (in)
dr = 6.00;              % Diameter at rear of transition (in)
Lt = 0;                 % Length of transition (in)
Xp = 0;                 % Distance from tip of nose to front of transition (in)
Cr = 10;                % Fin root chord (in)
Ct = 1.5;               % Fin tip chord (in)
S = 4.625;              % Fin semispan (in)
Lf = 5.58;              % Length of fin mid-chord line (in)
R = 2.00;               % Radius of body at aft end (in)
Xr = 7.125;             % Distance between fin root leading edge and fin tip leading edge parallel to body (in)
Xb = 119.5;             % Distance from nose tip to fin root chord leading edge (in)
N = 3;                  % Number of fins

%% Nose Cone Term

Xn = 0.666*Ln;          % For Ogive
Cnn = 2;

%% Conical Transition Terms
Cnt = 2*((dr/d)^2-(df/d)^2);
Xt = Xp+(Lt/3)*(1+((1-(df/dr))/(1-(df/dr)^2)));

%% Fin Terms

Cnf = (1+(R/(S+R)))*((4*N*(S/d)^2)/(1+(sqrt(1+(((2*Lf)/(Cr+Ct))^2)))));
Xf = Xb+((Xr*(Cr+2*Ct))/(3*(Cr+Ct)))+((1/6)*((Cr+Ct)-((Cr*Ct)/(Cr+Ct))));

%% Center of Pressure

Cnr = Cnn + Cnf;
X = ((Cnn*Xn)+(Cnf*Xf))/Cnr;               % Location from Nose Tip

fprintf('\n Center of Pressure Distance from Nose Cone         = %8.3f in', X)