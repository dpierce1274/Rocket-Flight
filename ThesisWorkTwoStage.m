%% Thesis Rocket
% The program Thesis Rocket is designed to calculate the apogee of a custom
% designed two-stage rocket with two O-motors. The motor thrust data used was from
% thrustcurve.com

%% Program Start
% Program start will clean the command window, workspace variables, and
% format Matlab output in a compact mode

clear, clc, format compact;

%% Given Rocket Characteristics
% Rocket Characteristics of the designed rocket

global Cd A rhos m1 O8000 g0 tb v pro array tend;

[datafile,path] = uigetfile({'*.csv'},'Select Data File');     % Selecting Flight Data
array = xlsread(datafile);      % Defining the array that contains the data

d = .1397;     % Body Tube Diameter of Rocket (m)
Cd = 0.72;     % Coefficient of Drag based on ogive nosecone
m1 = 46.35;    % Rocket Empty Mass (kg) (17.9 lbs)
O8000 = 14.062;% O8000 Rocket Engine Mass (kg) Empty Mass
pro = 18.61;   % O8000 Rocket Engine Propellant Mass (kg)

%% Constants Assumed

rhos = 1.225;       % Density of air in kg/m^3 at 20 deg. Centigrade
A = pi*(d*.5)^2;    % Cross-sectional Area of rocket (m^2)
g0 = 9.807;         % Gravity constant (m/s^2)

%% Boost Phase

tb = linspace(0,50,1000);       % Boost phase time
dt = tb(2)-tb(1);               % Time increment
v = zeros(1,size(tb,2));        % Empty set pre-allocation for velocity
a = zeros(1,size(tb,2));        % Empty set pre-allocation for acceleration
h = zeros(1,size(tb,2));        % Empty set pre-allocation for height
drag = zeros(1,size(tb,2));     % Empty set pre-allocation for drag
v(1) = 0;                       % Define v(1)
h(1) = 0;                       % Define h(1)
a(1) = 0;                       % Define a(1)
drag(1) = 0;                    % Define drag(1)
for k = 1:length(tb)-1
    a(k+1) = acceleration1(tb(k),v(k),h(k)); % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                % Euler's Method for ODE's to get height
    drag(k) = Drag(v(k),h(k));
end


%% Finding Thrust in Boost Phase 1

tpowered = round(tend/dt);  % Time in powered flight 

T = zeros(1,tpowered);      % Empty set pre-allocation for Thrust

for k = 1:tpowered
    T(k) = Thrust(tb(k));
end


%% Finding Mass in Boost Phase 1

m = zeros(1,tpowered);    % Empty set pre-allocation for Mass

for k = 1:tpowered
    m(k) = Mass1(tb(k));
end


%% Displaying Data
[z,y] = max(h);
disp(['Maximum height is ',num2str(z),'(m) at ',num2str(tb(y)),'(s)'])
[z,y] = max(v);
disp(['Maximum velocity is ',num2str(z),'(m/s) at ',num2str(tb(y)),'(s)'])
[z,y] = max(a);
disp(['Maximum acceleration is ',num2str(z),'(m/s^2) at ',num2str(tb(y)),'(s)'])
[z,y] = max(drag);
disp(['Maximum drag is ',num2str(z),'(m/s^2) at ',num2str(tb(y)),'(s)'])

%% Plotting the Data

figure(1)
plot(tb,v), xlabel('Time(s)'), ylabel('Velocity (m/s)'), title('Velocity vs. Time')
figure(2)
plot(tb,h), xlabel('Time(s)'), ylabel('Height (m)'), title('Height vs. Time')
figure(3)
plot(tb,a), xlabel('Time(s)'), ylabel('Acceleration (m/s^2)'), title('Acceleration vs. Time')
figure(4)
plot(tb(1:tpowered),T), xlabel('Time(s)'), ylabel('Thrust(N)'), title('Thrust in Boost Phase')
figure(5)
plot(tb(1:tpowered),m), xlabel('Time(s)'), ylabel('Mass(kg)'), title('Mass in Boost Phase')
figure(6)
plot(tb,drag), xlabel('Time(s)'), ylabel('Drag'), title('Drag vs. Time')