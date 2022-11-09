%% Parachute
% The Parachute script is used to calculate the descent rate of the various
% rocket components using a dual-deployment technique

%% Program Start
% Program start will clean the command window, workspace variables, and
% format Matlab output in a compact mode

clear, clc, format compact

%% Given Parachute Characteristics

d = 24;                % (Input) Parachute diam (in)
Cd = 1.55;             % (Input) Coefficient of Drag
mass = 25;             % (Input) Rocket Empty Mass (lbs)
Apogee = 250000;       % (Input) Max Altitude (ft)
main_Cd = 2;           % (Input) Main Parachute Cd
main_d = 90;           % (Input) Main Parachute diameter (in)

%% Constants Assumed

rhos = 1.225;       % Density of air in kg/m^3 at 20 deg. Centigrade
g0 = 9.807;         % Gravity constant (m/s^2)
Re = 6378;          % Radius of Earth (km)

%% Conversion Calculations

d = d*0.0254;                % Parachute diam (m)
mass = mass/2.2;             % Rocket Mass (kg)
h(1) = 0.3048*Apogee;        % Max Altitude (m)
A = pi*(d*.5)^2;             % Cross-sectional Area of paracute (m^2)
d_main = main_d*0.0254;      % Main Parachute diam(m)
A_main = pi*(d_main*.5)^2;   % Cross-sectional Area of main paracute (m^2)

%% Descent Calculations

rho(1) = rhos;
v(1) = 0;
gh(1) = g0;
dt = .1;
check = 0;

for k = 1:100000
    rho(k) = rhos*exp((-1/8000)*h(k));
    D(k) = .5*rho(k)*v(k)^2*Cd*A;
    dvdt(k) = (D(k)/mass)-g0;
    Vterm(k) = sqrt((2*mass*g0)/(rho(k)*A*Cd));
    if norm(v(k)) >= norm(Vterm)
        v(k) = Vterm(k);
        v(k+1) = Vterm(k);
        h(k+1) = h(k) + v(k)*dt;
    else
        v(k+1) = v(k) + dvdt(k)*dt;
        h(k+1) = h(k) + v(k)*dt;
    end
    if h(k) <= 457
        Cd = main_Cd;
        A = A_main;
        if check == 0
        fprintf('Descent Rate at Main Open = %2.2f m/s or %2.2f fps\n',v(k),v(k)*3.28084)
        check = 1;
        end
    end
    if h(k) <= 0
        fprintf('Descent Rate at Landing = %2.2f m/s or %2.2f fps\n',v(k),v(k)*3.28084)
        break
    end
end

t = 1:length(h);
fprintf('Total Elapsed Time = %.0f sec or %.2f min \n',t(end)/10,t(end)/600)
plot(t/10,h*3.28084/1000)
title('Altitude vs. Time'), xlabel('Time (s)'), ylabel('Altitude (ft)')
grid on
