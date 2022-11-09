%% Modified Barrowman Stability Equation

clear, clc, format compact


%% Input Parameters

R = 6/2;                    % Radius of nosecone at base (in)
L = 34.5;                   % Length of nosecone (in)
fidelity = 1000;            % Fidelity of calculation

%% Calculations

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



X = (L-(volume/(pi*y(end)^2)));

figure(1)
title('LD-Haack (Von Karman) Nosecone')
plot(x,y,x,neg_y,x(k),y(k),'or')
xlabel('Length (in)'), ylabel('Radius (in)')
axis('equal')

grid on

fprintf('LD-Haack Barrowman Parameter = %.3f \n',X/L)

    