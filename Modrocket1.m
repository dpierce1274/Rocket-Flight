%% EZI-65 (J90) Rocket
% The program EZI-65 (J90) Rocket is designed to calculate the apogee of an 
% LOC EZI-65 rocket with a J90 motor. The motor thrust data used was from
%

%% Program Start
% Program start will clean the command window, workspace variables, and
% format Matlab output in a compact mode

clear, clc, format compact;

%% Given Rocket Characteristics
% Rocket Characteristics that are given in the problem

global Cd A rhos m1 J90 g0 tb v pro array;

[datafile,path] = uigetfile({'*.csv'},'Select Data File');     % Selecting Flight Data
array = xlsread(datafile);      % Defining the array that contains the data

d = .1016;    % Body Tube Diameter of Rocket (m)
Cd = 0.36;    % Coefficient of Drag based on ogive nosecone
m1 = 1.244;   % Rocket Empty Mass (kg)
J90 = 0.443;  % J90 Rocket Engine Mass (kg) Empty Mass
pro = .391;   % J90 Rocket Engine Propellant Mass (kg)

%% Constants Assumed

rhos = 1.225;       % Density of air in kg/m^3 at 20 deg. Centigrade
A = pi*(d*.5)^2;    % Cross-sectional Area of rocket (m^2)
g0 = 9.807;         % Gravity constant (m/s^2)

%% Boost Phase

tb = linspace(0,18,1000);       % Boost phase time
dt = tb(2)-tb(1);               % Time increment
v = zeros(1,size(tb,2));        % Empty set pre-allocation for velocity
a = zeros(1,size(tb,2));        % Empty set pre-allocation for acceleration
h = zeros(1,size(tb,2));        % Empty set pre-allocation for height
v(1) = 0;                       % Define v(1)
h(1) = 0;                       % Define h(1)
a(1) = 0;                       % Define a(1)
for k = 1:length(tb)-1
    a(k+1) = acceleration(tb(k),v(k),h(k)); % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                % Euler's Method for ODE's to get height
end


%% Finding Thrust in Boost Phase

T = zeros(1,500);    % Empty set pre-allocation for Thrust

for k = 1:500
    T(k) = Thrust(tb(k));
end


%% Finding Mass in Boost Phase

m = zeros(1,500);    % Empty set pre-allocation for Mass

for k = 1:500
    m(k) = Mass(tb(k));
end


%% Displaying Data
[z,y] = max(h);
disp(['Maximum height is ',num2str(z),'(m) at ',num2str(tb(y)),'(s)'])
[z,y] = max(v);
disp(['Maximum velocity is ',num2str(z),'(m/s) at ',num2str(tb(y)),'(s)'])
[z,y] = max(a);
disp(['Maximum acceleration is ',num2str(z),'(m/s^2) at ',num2str(tb(y)),'(s)'])

%% Plotting the Data

figure(1)
plot(tb,v), xlabel('Time(s)'), ylabel('Velocity (m/s)'), title('Velocity vs. Time')
figure(2)
plot(tb,h), xlabel('Time(s)'), ylabel('Height (m)'), title('Height vs. Time')
figure(3)
plot(tb,a), xlabel('Time(s)'), ylabel('Acceleration (m/s^2)'), title('Acceleration vs. Time')
figure(4)
plot(tb(1:500),T), xlabel('Time(s)'), ylabel('Thrust(N)'), title('Thrust in Boost Phase')
figure(5)
plot(tb(1:500),m), xlabel('Time(s)'), ylabel('Mass(kg)'), title('Mass in Boost Phase')




