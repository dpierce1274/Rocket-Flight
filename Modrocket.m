%% Astro Modrocket
% The program Modrocket is designed to calculate the apogee of an Estes Model
% Rocket

%% Program Start
% Program start will clean the command window, workspace variables, and
% format Matlab output in a compact mode

clear, clc, format compact;

%% Given Rocket Characteristics
% Rocket Characteristics that are given in the problem%% Select and Import Data File

[datafile,path] = uigetfile({'*.csv'},'Select Data File');     % Selecting Flight Data
array = xlsread(datafile);      % Defining the array that contains the data
tb = transpose((array(1:end,1)));      % Burn Time Data from text file (s)
Tdata = transpose((array(1:end,2)));  % Thrust Data from text file (N)
Thrust = interp1(tb,Tdata,'linear')

%% Assigning Data

d = .1016;             % (Input) Body Tube Diameter of Rocket (m)
Cd = 0.75;             % (Input) Coefficient of Drag
m1 = 1.244;            % (Input) Rocket Empty Mass (kg)
J90 = 0.443;           % J90 Rocket Engine Mass (kg) Empty Mass
pro = .391;            % J90 Rocket Engine Propellant Mass (kg)
mint = m1 + J90 + pro; % Initial Mass of Loaded Rocket (kg)
%% Constants Assumed

rhos = 1.225;       % Density of air in kg/m^3 at 20 deg. Centigrade
A = pi*(d*.5)^2;    % Cross-sectional Area of rocket (m^2)
g0 = 9.807;         % Gravity constant (m/s^2)
Re = 6378;          % Radius of Earth (km)

%% Boost Phase
rho(1) = rhos;
h(1) = 0;
v(1) = 0;
gh(1) = g0;
dmdt = pro/tb(end);
dt(1) = 1;
Mass(1) = mint;
time = tb;

for k = 1:100
    dt(k+1) = tb(k+1)-tb(k);
    rho(k) = rhos*exp((-1/8000)*h(k));
    Drag(k+1) = .5*rho(k)*v(k)^2*Cd*A;
    Mass(k+1) = -dmdt*tb(k)+mint;
    gh(k+1) = g0*(Re/(Re+h(k)))^2;
    dvdt(k) = (Thrust(k)-Drag(k))/Mass(k)-gh(k); % Acceleration function to find dvdt
    v(k+1) = v(k) + dvdt(k)*dt(k);                 % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt(k);
end

%% Coast Phase

timecoast = v(length(100))/g0;
dtcoast = timecoast/500;
endcoast = 500+length(100);

for k = 101:endcoast
    time(k+1) = time(k)+dtcoast;
    rho(k) = rhos*exp((-1/8000)*h(k));
    Drag(k+1) = .5*rho(k)*v(k)^2*Cd*A;
    Mass(k+1) = mint;
    gh(k+1) = g0*(Re/(Re+h(k)))^2;
    dvdt(k) = -Drag(k)/Mass(k)-gh(k);
    v(k+1) = v(k) + dvdt(k)*dtcoast;                 % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dtcoast;
end


%% Finding Height
% Height is determined using the cumtrapz function to intergrate velocity


h = cumtrapz(time,v);   % Integration of velocity to get height (m)

%% Displaying Data
[z,y] = max(h);
disp(['Maximum height is ',num2str(z),'(m) at ',num2str(time(y)),'(s)'])
[z,y] = max(v);
disp(['Maximum velocity is ',num2str(z),'(m/s) at ',num2str(time(y)),'(s)'])
%[z,y] = max(a);
%disp(['Maximum acceleration is ',num2str(z),'(m/s^2) at ',num2str(time(y)),'(s)'])

%% Plotting the Data



figure(1)
plot(time,v), xlabel('Time(s)'), ylabel('Velocity (m/s)'), title('Velocity vs. Time')
figure(2)
plot(time,h), xlabel('Time(s)'), ylabel('Height (m)'), title('Height vs. Time')
figure(3)
%plot(time,dvdt), xlabel('Time(s)'), ylabel('Acceleration (m/s^2)'), title('Acceleration vs. Time')
%figure(4)
%plot(tb,T), xlabel('Time(s)'), ylabel('Thrust(N)'), title('Thrust in Boost Phase')
%figure(5)
%plot(tb,m), xlabel('Time(s)'), ylabel('Mass(kg)'), title('Mass in Boost Phase')


