%% Thesis Rocket
% The program ThesisRocketSustainerOnly is designed to calculate key flight 
% performance parameters of the NPS Class-3 High Power Rocket 
% (sustainer-only configuration). Motor thrust data provided by 
% thrustcurve.com for the O3400 sustainer motor.

%The example output for this specific script shows the estimated key flight 
%performance parameters for the single-stage configuration using the
%redesigned sustainer (flight 6).

%Source Author: Dillon Pierce
%Name of File: ThesisRocketSustainerOnly.m
%File Location: ssagcommon$\High Power Rocket\MATLAB
%Date Last Modified: 23 May 2019

%Inputs: Various rocket parameters
%Outputs: Rocket key flight performance parameters

%% Program Start
% Program start will clean the command window, workspace variables, and
% format Matlab output in a compact mode

clear, clc, format compact;

%% Given Rocket Characteristics and Other Inputs
% Defined Rocket Characteristics

global rhos g0

[sustainer_datafile,~] = uigetfile({'*.csv'},'Select Sustainer Thrust Data File');  % Select Sustainer Thrust Data
sustainer_thrust_array = xlsread(sustainer_datafile);      % Defining the array that contains the data

d_sustainer = .1016;            % Body Tube Diameter of Sustainer (m)
Cd = 0.7;                       % Coefficient of Drag based on Von Karman nosecone (.516 from test data)
m_sustainer = 12;               % Sustainer Empty Mass (kg) (25 lbs)
m_O3400_tot = 16.84;            % Sustainer Motor Total Mass (kg) (37 lbs)
m_O3400_prop = 10.93;           % Sustainer Motor Propellant Mass (kg)(24 lbs)

A_sust = pi*(d_sustainer*.5)^2; % Cross-sectional Area of sustainer (m^2)

fidelity = 1000;    % Fidelity increment for calculations
launch_alt = 626;   % Launch site altitude (m)

%% Constants Assumed

rhos = 1.225;       % Density of air in kg/m^3 at 20 deg. Centigrade
g0 = 9.807;         % Gravity constant (m/s^2)



%% Preset Allocations

v(1) = 0;                          % Define v(1)
h(1) = launch_alt;                 % Define h(1) (Launch Site Altitude (m))
a(1) = 0;                          % Define a(1)
rho(1) = rhos*exp((-1/7800)*h(1)); % Define rho(1)


%% Sustainer Thrust Phase

% Calculate Thrust Array from Excel Spreadsheet
t_sust = linspace(sustainer_thrust_array(1,1),sustainer_thrust_array(end,1),fidelity);   % Boost phase time
tdata = sustainer_thrust_array(1:end,1);                        % Time Data from excel spreadsheet (s)
Tdata = sustainer_thrust_array(1:end,2);                        % Thrust Data from excel spreadsheet (N)
T_sust = interp1(tdata,Tdata,t_sust,'linear');    % Interpolation to get thrust (N) at each time increment
t_boost_sust = t_sust;
t = t_boost_sust;

% Calculate Mass Array of rocket during sustainer thrust phase
mint = m_sustainer + m_O3400_tot;                % Initial Mass (kg)  (Empty stage + motor (kg))
mfin = m_sustainer + m_O3400_tot - m_O3400_prop; % Final Mass (kg) (Empty stage + sustainer casing (kg))
dmdt = (mint-mfin)/(t_sust(end));  % Linear Mass change over burnout time
mass_sustboost = -t_sust.*dmdt+mint;
mass = mass_sustboost;

dt = t_sust(2)-t_sust(1);     % Time increment

for k = 1:length(t)-1
    rho(k+1) = Density(h(k));                    % Atmoshperic Density
    D(k+1) = Drag(v(k),h(k),A_sust,Cd);         % Function for drag (N)
    a(k+1) = acceleration(T_sust(k),D(k),mass_sustboost(k)); % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                     % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                     % Euler's Method for ODE's to get height
end

%% Sustainer Coast phase

apogee = 50;
t_coast2 = linspace(t(end),apogee,fidelity); % Coast Time
dt = t_coast2(2)-t_coast2(1);
t_coast2 = t_coast2+dt;                      % Coast Time

for k = length(t)-1:length(t)+length(t_coast2)-1
    mass(k+1) = mfin;    
    rho(k+1) = Density(h(k));                    % Atmoshperic Density
    D(k+1) = Drag(v(k),h(k),A_sust,Cd);          % Function for drag (N)
    a(k+1) = acceleration(0,D(k),mfin);          % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                     % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                     % Euler's Method for ODE's to get height
end

t = horzcat(t,t_coast2);

%% Dynamic Pressure

for k = 1:length(v)
    q(k) = .5*rho(k)*v(k)^2;        % Aerodynamic Pressure (Pa)
end


%% Displaying Data
[z,y] = max(h);
fprintf('Maximum height is %0.f (m) or %0.f (ft) at %.0f (s)\n',z,z*3.281,t(y))
[z,y] = max(v);
fprintf('Maximum velocity is %0.f (m/s) or %0.f (ft/s) at %.0f (s)\n',z,z*3.281,t(y))
[z,y] = max(a);
fprintf('Maximum acceleration is %0.f (m/s^2) or %0.f (ft/s^2) at %.0f (s)\n',z,z*3.281,t(y))
[maxq,index_maxq] = max(q);
fprintf('Maximum Dynamic Pressure: %0.f (Pa) at %.0f (s) and %0.f (m-AGL)\n',maxq,t(index_maxq),h(index_maxq))
%% Plotting the Data

figure(1)
plot(t,a), xlabel('Time(s)'), ylabel('Acceleration (m/s^2)'), title('Acceleration vs. Time')
grid on
figure(2)
plot(t,v), xlabel('Time(s)'), ylabel('Velocity (m/s)'), title('Velocity vs. Time')
grid on
figure(3)
plot(t,h), xlabel('Time(s)'), ylabel('Height (m)'), title('Height vs. Time')
grid on
figure(4)
plot(t,mass), xlabel('Time(s)'), ylabel('Mass (kg)'), title('Rocket Mass vs. Time')
grid on
figure(5)
plot(t_sust,T_sust), xlabel('Time(s)'), ylabel('Thrust(N)'), title('Thrust in Boost Phase')
grid on
figure(6)
plot(t,rho), xlabel('Time(s)'), ylabel('Air Density (kg/m^3)'), title('Air Density vs. Time')
grid on
figure(7)
subplot(2,1,1)
plot(t,q,t(index_maxq),maxq,'or'), xlabel('Time(s)'), ylabel('Aerodynamic Pressure (Pa)'), title('Dynamic Pressure vs. Time')
legend('Dynamic Pressure (Pa)','Max Q')
grid on
subplot(2,1,2)
plot(h,q,h(index_maxq),maxq,'or'), xlabel('Altitude (m)'), ylabel('Aerodynamic Pressure (Pa)'), title('Dynamic Pressure vs. Altitude')
legend('Dynamic Pressure (Pa)','Max Q')
grid on
%% Functions

function [ D ] = Drag( v,h,A,Cd )
%The Drag function calculates the drag on an object at a certain velocity

D = .5*Density(h)*v^2*Cd*A;   % Function for drag (N)


end

function [rho] = Density(h)
% Function for changing atmospheric density
global rhos;

rho = rhos*exp((-1/7800)*h); 

end
    
function [ dvdt ] = acceleration(T,D,M)
% The acceleration function is used to calculate the acceleration of a
% rocket at a given time and velocity.

global g0;
 
dvdt = (T - D)/M - g0;



end
