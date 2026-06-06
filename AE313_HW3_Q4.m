% Description: This code will use the ODE45 method to solve the second
% order diferential equation of a two-body motion about earth and fidn
% orbital period, trajectory and initial position.
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE313 - Section 01DB
% Date: 09/30/2025
% Worked With: N/a

% Iniciation
clear;
clc;
close ALL;

%set inputs
mu = 398600; %Earth's Gravitational parameter
r0 = [7000; 0; 0]; %given postion vector
v0 = [0.2; 7; 4];   %given velocity vector
Y0 = [r0; v0];   %initial condition for ODE45 

%Orbital period calculator for comparison
r0_mag = norm(r0);
V0_mag = norm(v0);
epsilon = (V0_mag^2 / 2) - (mu / r0_mag);
a = -mu / (2 * epsilon);
T = 2 * pi * sqrt(a^3 / mu); 

%create two-body propagator fuction iterator
function dYdt = twoBodyPropagator(~, Y, mu)
    r_vec = Y(1:3);
    r = norm(r_vec);
    acceleration = -mu / (r^3) * r_vec;
    dYdt = [Y(4:6); acceleration];
end

%Propagation of Orbit using ode45
tspan = [0, T];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); %define accuracy
[t, Y] = ode45(@(t, Y) twoBodyPropagator(t, Y, mu), tspan, Y0, options);

%Plot of trajectory
figure('Color', 'w');
plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'LineWidth', 1.5);
hold on;

% Mark the initial position
plot3(Y(1, 1), Y(1, 2), Y(1, 3), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(Y(1, 1), Y(1, 2), Y(1, 3), '  r_0', 'FontSize', 10);
title('3D plot of orbit with intial position');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
grid on;
axis equal; 
view(30, 30);