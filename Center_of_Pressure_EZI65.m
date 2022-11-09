%% Computing Center of Pressure (EZI-65)

clear, clc, format compact

%% Inputs

Ln = 12.75;             % Length of nosecone (in)
d = 4.00;               % Diameter at base of nose (in)
df = 4.00;              % Diameter at front of transition (in)
dr = 4.00;              % Diameter at rear of transittion (in)
Lt = 0;                 % Length of transition (in)
Xp = 0;                 % Distance from tip of nose to front of transition (in)
Cr = 5.50;              % Fin root chord (in)
Ct = 3.5;               % Fin tip chord (in)
S = 4.75;               % Fin semispan (in)
Lf = 4.75;              % Length of fin mid-chord line (in)
R = 2.00;               % Radius of body at aft end (in)
Xr = 1.00;              % Distance between fin root leading edge and fin tip leading edge parallel to body (in)
Xb = 52;                % Distance from nose tip to fin root chord leading edge (in)
N = 3;                  % Number of fins

%% Nose Cone Term

Xn = 0.466*Ln;          % For Ogive
Cnn = 2;


%% Fin Terms

Cnf = (1+(R/(S+R)))*((4*N*(S/d)^2)/(1+(sqrt(1+(((2*Lf)/(Cr+Ct))^2)))));
Xf = Xb+((Xr*(Cr+2*Ct))/(3*(Cr+Ct)))+((1/6)*((Cr+Ct)-((Cr*Ct)/(Cr+Ct))));

%% Center of Pressure

Cnr = Cnn + Cnf;
X = ((Cnn*Xn)+(Cnf*Xf))/Cnr;               % Location from Nose Tip

fprintf('\n Center of Pressure Distance from Nose Cone         = %8.3f in', X)