%% Thesis Rocket
% The program Thesis Rocket is designed to calculate key flight performance 
% parameters of the NPS Class-3 High Power Rocket. Motor thrust data
% provided by thrustcurve.com for the O3400 sustainer motor and Robert
% DeHate of Animal Motor Works for the Q15782 booster motor.

%% Program Start
% Program start will clean the command window, workspace variables, and
% format Matlab output in a compact mode

clear, clc, format compact;

%% Given Rocket Characteristics and Other Inputs
% Defined Rocket Characteristics

global rhos g0

[booster_datafile,~] = uigetfile({'*.csv'},'Select Booster Thrust Data File');     % Select Booster Thrust Data
booster_thrust_array = xlsread(booster_datafile);      % Defining the array that contains the data
[sustainer_datafile,~] = uigetfile({'*.csv'},'Select Sustainer Thrust Data File');  % Select Sustainer Thrust Data
sustainer_thrust_array = xlsread(sustainer_datafile);      % Defining the array that contains the data

d_booster = .1524;              % Body Tube Diameter of Booster (m)
d_sustainer = .1016;            % Body Tube Diameter of Sustainer (m)
Cd = 0.7;                       % Coefficient of Drag based on Von Karman nosecone (.516 from test data)
m_booster = 4.6;                % Booster Empty Mass (kg) (10 lbs)
m_sustainer = 12;               % Sustainer Empty Mass (kg) (12 lbs)
m_tot = m_booster+m_sustainer;  % Total Rocket mass (kg)
m_Q15782_prop = 40.75;          % Booster Motor Propellant Mass (kg) (90 lbs)
m_Q15782_tot = 49.84;           % Booster Motor Total Mass (kg) (110 lbs)
m_O3400_tot = 16.84;            % Sustainer Motor Total Mass (kg) (37 lbs)
m_O3400_prop = 10.93;           % Sustainer Motor Propellant Mass (kg)(24 lbs)

A_boost = pi*(d_booster*.5)^2;  % Cross-sectional Area of booster (m^2)
A_sust = pi*(d_sustainer*.5)^2; % Cross-sectional Area of sustainer (m^2)

fidelity = 1000;    % Fidelity increment for calculations
launch_alt = 626;   % Launch site altitude (m)
stage_delay = 17;   % Sustainer Stage ignition delay from launch (s)
sep_time = 17;      % Estimated time of stage seperation (s)

%% Constants Assumed

rhos = 1.225;       % Density of air in kg/m^3 at 20 deg. Centigrade
g0 = 9.807;         % Gravity constant (m/s^2)


%% Booster Thrust Phase

% Calculate Thrust Array from Excel Spreadsheet
t_boost = linspace(booster_thrust_array(1,1),booster_thrust_array(end,1),fidelity);   % Boost phase time
tdata = booster_thrust_array(1:end,1);                        % Time Data from excel spreadsheet (s)
Tdata = booster_thrust_array(1:end,2);                        % Thrust Data from excel spreadsheet (N)
T_boost = interp1(tdata,Tdata,t_boost,'linear');    % Interpolation to get thrust (N) at each time increment


% Calculate Mass Array of rocket during booster thrust phase
mint = m_tot+m_Q15782_tot+m_O3400_tot;               % Initial Mass (kg)  (Empty Rocket + Motors)
mfin = m_tot+m_O3400_tot+m_Q15782_tot-m_Q15782_prop; % Final Mass (kg) (Empty Rocket + sustainer motor - empty booster)
dmdt = (mint-mfin)/t_boost(end);  % Linear Mass change over burnout time
mass_boostphase = -t_boost.*dmdt+mint;
mass = mass_boostphase;


% Preset Allocations
dt = t_boost(2)-t_boost(1);        % Time increment
v = zeros(1,length(t_boost));      % Empty set pre-allocation for velocity
a = zeros(1,length(t_boost));      % Empty set pre-allocation for acceleration
h = zeros(1,length(t_boost));      % Empty set pre-allocation for height
rho = zeros(1,length(t_boost));    % Empty set pre-allocation for Density
D = zeros(1,length(t_boost));      % Empty set pre-allocation for Drag
v(1) = 0;                          % Define v(1)
h(1) = launch_alt;                 % Define h(1) (Launch Site Altitude (m))
a(1) = 0;                          % Define a(1)
rho(1) = rhos*exp((-1/7800)*h(1)); % Define rho(1)


for k = 1:length(t_boost)-1
    rho(k+1) = Density(h(k));                    % Atmoshperic Density
    D(k+1) = Drag(v(k),h(k),A_boost,Cd);         % Function for drag (N)
    a(k+1) = acceleration(T_boost(k),D(k),mass_boostphase(k)); % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                     % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                     % Euler's Method for ODE's to get height
end


%% Coast phase
t_coast1 = linspace(t_boost(end),stage_delay,fidelity); % Coast Time
t = horzcat(t_boost,t_coast1);

for k = length(t_boost):length(t_boost)+length(t_coast1)-1
    if t_coast1(k-length(t_boost)+1) > sep_time
        A = A_sust;
        Cd = Cd;
        m = m_sustainer + m_O3400_tot;
    else
        A = A_boost;
        Cd = Cd;
        m = mfin;
    end 
    mass(k+1) = m;    
    rho(k+1) = Density(h(k));                    % Atmoshperic Density
    D(k+1) = Drag(v(k),h(k),A,Cd);               % Function for drag (N)
    a(k+1) = acceleration(0,D(k),m);             % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                     % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                     % Euler's Method for ODE's to get height
end

%% Sustainer Thrust Phase

% Calculate Thrust Array from Excel Spreadsheet
t_sust = linspace(sustainer_thrust_array(1,1),sustainer_thrust_array(end,1),fidelity);   % Boost phase time
tdata = sustainer_thrust_array(1:end,1);                        % Time Data from excel spreadsheet (s)
Tdata = sustainer_thrust_array(1:end,2);                        % Thrust Data from excel spreadsheet (N)
T_sust = interp1(tdata,Tdata,t_sust,'linear');    % Interpolation to get thrust (N) at each time increment
t_boost_sust = t_sust+t_coast1(end);
t = horzcat(t,t_boost_sust);

% Calculate Mass Array of rocket during sustainer thrust phase
mint = m_sustainer + m_O3400_tot;                % Initial Mass (kg)  (Empty stage + motor (kg))
mfin = m_sustainer + m_O3400_tot - m_O3400_prop; % Final Mass (kg) (Empty stage + sustainer casing (kg))
dmdt = (mint-mfin)/(t_sust(end));  % Linear Mass change over burnout time
mass_sustboost = -t_sust.*dmdt+mint;
mass = horzcat(mass,mass_sustboost);

dt = t_sust(2)-t_sust(1);     % Time increment

for k = length(t_boost)+length(t_coast1):length(t)-1
    rho(k+1) = Density(h(k));                    % Atmoshperic Density
    D(k+1) = Drag(v(k),h(k),A_sust,Cd);         % Function for drag (N)
    a(k+1) = acceleration(T_sust(k-(fidelity*2)+1),D(k),mass_sustboost(k-(fidelity*2-1))); % Acceleration function to find dvdt
    v(k+1) = v(k) + a(k)*dt;                     % Euler's Method for ODE's to get velocity
    h(k+1) = h(k) + v(k)*dt;                     % Euler's Method for ODE's to get height
end

%% Sustainer Coast phase
apogee = 150;
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
plot(t,v)


%% Displaying Data
[z,y] = max(h);
fprintf('Maximum height is %0.f (m) or %0.f (ft) at %.0f (s)\n',z,z*3.281,t(y))
[z,y] = max(v);
fprintf('Maximum velocity is %0.f (m/s) or %0.f (ft/s) at %.0f (s)\n',z,z*3.281,t(y))
[z,y] = max(a);
fprintf('Maximum acceleration is %0.f (m/s^2) or %0.f (ft/s^2) at %.0f (s)\n',z,z*3.281,t(y))
fprintf('Velocity at sustainer stage ignition: %0.f (m/s) or %0.f (ft/s) at %.0f (s)\n',v(2000),v(2000)*3.281,t(y))
%% Plotting the Data

figure(1)
plot(t,a), xlabel('Time(s)'), ylabel('Acceleration (m/s^2)'), title('Acceleration vs. Time')
figure(2)
plot(t,v), xlabel('Time(s)'), ylabel('Velocity (m/s)'), title('Velocity vs. Time')
figure(3)
plot(t,h), xlabel('Time(s)'), ylabel('Height (m)'), title('Height vs. Time')
figure(4)
plot(t_boost,T_boost), xlabel('Time(s)'), ylabel('Thrust(N)'), title('Thrust in Booster Thrust Phase')
figure(5)
plot(t_boost_sust,T_sust), xlabel('Time(s)'), ylabel('Thrust(N)'), title('Thrust in Sustainer Thrust Phase')
figure(6)
plot(t,mass), xlabel('Time(s)'), ylabel('Mass (kg)'), title('Rocket Mass vs. Time')
figure(7)
plot(t,rho), xlabel('Time(s)'), ylabel('Air Density (kg/m^3)'), title('Air Density vs. Time')


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

