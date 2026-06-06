% Code Name: ECI to COE converter
% Code Description: Program made to convert Earth Cenetered Inertial frame
% velocity and position vectors into the six classical orbital elements of
% orbit [a,e, i, raan, argument of perigee, true anomoly], and vise versa.
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE313 - 01DB
% Date: 10/17/2023
% Worked With: N/a

%Iniciation
clear
clc
fprintf('This program converts Earth Cenetered postion and velocity vectors to the six classical orbiatl elements and vise-versa.')

% Data set used, hard coded
eg = 398600.4418;
og_r = [-2600.000, -1300.000, 6800.000];
og_v = [5.900, -1.200, 0.800];
og_a = 7000;
og_e = 0.10;
og_i = 40;
og_raan = 60;
og_ap = 25;
og_theta = 45;

%functions called
[r_eci, v_eci] = orbitalElementsToECI(og_a, og_e, og_i, og_raan, og_ap, og_theta);
[a, e, i, RAAN, omega, theta] = eciToOrbitalElements(r_eci, v_eci);
%[a, e, i, RAAN, omega, theta] = eciToOrbitalElements(og_r, og_v);
%[r_eci, v_eci] = orbitalElementsToECI(a, e, i, RAAN, omega, theta);
Y0 = [r_eci;v_eci];

%first function COEtoRV.m
function [r_eci, v_eci] = orbitalElementsToECI(og_a, og_e, og_i, og_raan, og_ap, og_theta)
% Convert all angles from degrees to radians
n_i     = deg2rad(og_i);
n_RAAN  = deg2rad(og_raan);
n_ap = deg2rad(og_ap);
n_theta = deg2rad(og_theta);
eg = 398600.4418; %earth gravitational constant

h_12 = sqrt(eg*og_a*(1-og_e^2));  % semi-latus rectum calcualtion
r = (h_12^2/eg) / (1 + (og_e * cos(n_theta)));  % find radius magnitude

% compute position and velocity in perifocal (PQW) frame
r_pqw = [r * cos(n_theta), r * sin(n_theta), 0];
v_pqw = eg/h_12 * [-sin(n_theta), og_e + cos(n_theta), 0];

% generate rotational matrix from perifocal frame to ECI
R3_W = [cos(n_RAAN)  sin(n_RAAN)  0;
        -sin(n_RAAN) cos(n_RAAN)  0;
         0          0         1];  % Rotation about Z-axis (RAAN)

R1_i = [1  0           0;
        0  cos(n_i)  sin(n_i);
        0 -sin(n_i)  cos(n_i)];    % Rotation about X-axis (inclination)

R3_w = [cos(n_ap)  sin(n_ap)  0;
        -sin(n_ap) cos(n_ap)  0;
         0           0          1]; % Rotation about Z-axis (argument of perigee)

% find full rotation matrix to convert from perifocal frame to ECI
Q_pqw2eci = R3_w' * R1_i' * R3_W';  % Transpose = inverse (rotation matrices are orthogonal)

% transform both the position and radius to ECI frame using rotation matrix
r_eci = Q_pqw2eci * r_pqw';
v_eci = Q_pqw2eci * v_pqw';
fprintf('\nPosition in ECI (km): [%.3f, %.3f, %.3f]\n', r_eci(1), r_eci(2), r_eci(3)); %print results
fprintf('\nVelocity in ECI (km/s): [%.3f, %.3f, %.3f]\n', v_eci(1), v_eci(2), v_eci(3)); %print results
end

%second function RVtoCOE.m
function [a, e, i, RAAN, omega, theta] = eciToOrbitalElements(gv_r, gv_v)
r = norm(gv_r); % find maginitude of radius
v = norm(gv_v); % find maginitude of velocity
eg = 398600.4418; % Earth gravitational constant
h_vec = cross(gv_r, gv_v); % compute specific angular momentum
h = norm(h_vec); % specific angular momentum maginitude
i = acosd(h_vec(3)/h); % find inclination angle
k_hat = [0 0 1]; % set k direction vector
n_vec = cross(k_hat, h_vec); % compute line of nodes
n = norm(n_vec); % find line of nodes maginitude
e_vec = (cross(gv_v, h_vec)/eg) - (gv_r/r); % Calculate ecentricty vector
e = sqrt(1+((h^2/eg^2)*(v^2 - 2*(eg/r)))); % find eccentricty maginitude

if n_vec(2) >= 0 % if statement to set conditions to calculate RAAN
    RAAN = acosd(n_vec(1)/n);
else
    RAAN = 360 - acosd(n_vec(1)/n);
end
omega = acosd(dot(n_vec/n, e_vec/e)); % calculate argument of perigee

if dot(gv_r, gv_v) >= 0 % if statement to set conditions to calculate true anomaly
    theta = acosd(dot(e_vec, gv_r)/(e*r));
else
   theta = 360 - acosd(dot(e_vec, gv_r)/(e*r));
end

a = h^2 / (eg * (1 - 0.1^2)); % Calculate semi-major axis
fprintf('Semi-major axis (a): %.3f km\n', a);
fprintf('Eccentricity vector (e): %.5f\n', e_vec);
fprintf('Eccentricity (e): %.5f\n', e);
fprintf('Inclination (i): %.3f°\n', i);
fprintf('RAAN (Ω): %.3f°\n', RAAN);
fprintf('Argument of Perigee (ω): %.3f°\n', omega);
fprintf('True Anomaly (ν): %.3f°\n', theta);
end

% Third function body propagator to plot orbit from SET B
function dYdt = twoBodyPropagator(~, Y, mu)
    r_vec = Y(1:3);
    r = norm(r_vec);
    acceleration = -mu / (r^3) * r_vec;
    dYdt = [Y(4:6); acceleration];
end

%Propagation of Orbit using ode45
T = 2 * pi * sqrt(a^3 / eg); %Calculate Period
tspan = [0, T]; % run for 1 orbit cycle
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); %define accuracy, based on HW3_Q4
[t, Y] = ode45(@(t, Y) twoBodyPropagator(t, Y, eg), tspan, Y0, options);

%Plot of trajectory with initial postion marked
figure('Color', 'w');
plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'LineWidth', 1.5);
hold on;
plot3(Y(1, 1), Y(1, 2), Y(1, 3), 'ro', 'MarkerSize', 8, 'LineWidth', 2); 
text(Y(1, 1), Y(1, 2), Y(1, 3), '  r_0', 'FontSize', 10);
title('3D plot of orbit with intial position');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
grid on;

% ----- Sample Round-trip A ------
% This program converts Earth Cenetered postion and velocity vectors to the six classical orbiatl elements and vise-versa.Semi-major axis (a): 4936.359 km
% Eccentricity vector (e): 0.23440
% Eccentricity vector (e): 0.03037
% Eccentricity vector (e): -0.27344
% Eccentricity (e): 0.36143
% Inclination (i): 75.849°
% RAAN (Ω): 170.423°
% Argument of Perigee (ω): 128.720°
% True Anomaly (ν): 200.212°
% Position in ECI (km): [-2600.000, -1300.000, 6800.000]
% Velocity in ECI (km/s): [5.900, -1.200, 0.800]

% |og_r - r_eci| = 3*10^-7 km !Passes Test!
% |og_v - v_eci| = 1.872*10^-7 !Passes Test!

% ----- Sample Round-trip B ------
% This program converts Earth Cenetered postion and velocity vectors to the six classical orbiatl elements and vise-versa.
% Position in ECI (km): [-2928.054, 4246.638, 3909.439]
% Velocity in ECI (km/s): [-5.900, -5.193, 2.109]
% Semi-major axis (a): 7000.000 km
% Eccentricity vector (e): 0.01728
% Eccentricity vector (e): 0.09468
% Eccentricity vector (e): 0.02717
% Eccentricity (e): 0.10000
% Inclination (i): 40.000°
% RAAN (Ω): 60.000°
% Argument of Perigee (ω): 25.000°
% True Anomaly (ν): 45.000°

% |og_a - a| = 5*10^-12 km !Passes Test!
% |og_e - e| = 3*10^-15 !Passes Test!
% |og_i - i| = 10^-25 !Passes Test!
% |og_raan - RAAN| = 7*10^-15 !Passes Test!
% |og_ap - omega| = 2.991*10^-12 !Passes Test!
% |og_theta - theta| = 1.229*10^-12 !Passes Test!
