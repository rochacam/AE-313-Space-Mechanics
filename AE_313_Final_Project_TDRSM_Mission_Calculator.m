% Code Name: Final Project TDRS_M mission calculator
% Code Description: Program created to perfurm the orbital mechanics pahses
% of the TDRS_M mission computing the total delta velocity required to
% place satellite into specific communication GEO orbit from inital launch and placement in LEO orbit
% Authors: Matheus Rocha Carlos and
% Emails: ROCHACAM@my.erau.edu
% Class: AE313 - 01DB
% Dates worked: 11/10/2025 - 11/15/2025
% Outside Work: MATLAB Data Base and https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html

%Iniciation
clear; close all; clc;
% Constants
mu = 398600.4418;        % Earth's gravitational parameter [km^3/s^2]
RE = 6378.137;           % Earth's mean radius [km]
deg = pi/180;            % Degree to radians converter

%% Orbit Propagator (Mission Phase 1)
% Given parking orbit classical orbital elements
a = RE + 300;            % semi-major axis [km]
e = 0.001;               % eccentricity
i = 28.5 * deg;          % inclination [rad]
RAAN = 40 * deg;         % right ascension of ascending node [rad]
omega = 0 * deg;         % argument of perigee [rad]
nu0 = 0 * deg;           % initial true anomaly [rad]

% Algorithem one - convert Classical orbital elements to position and
% velocity vectors in ECI
[r0, v0] = coe2rv(a, e, i, RAAN, omega, nu0, mu);

% Algorithem two - propagate one full revolution using two-body motion
T = 2*pi*sqrt(a^3/mu); % compute orbital period (1)
y0 = [r0; v0]; % get initial state vector from algorithem one

% ODE45 options and integration
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
nplot = 100000; tspan = linspace(0, T, nplot); % timesetp
[tt, yy] = ode45(@(t,y) twoBodyEOM(t,y,mu), tspan, y0, opts);

%get r and v vectors and respective maginitude history during propagation
r_hist = yy(:,1:3);
v_hist = yy(:,4:6);
r_mag = sqrt(sum(r_hist.^2,2));
v_mag = sqrt(sum(v_hist.^2,2));

% Compute COE changes vs time from propagation
n_steps = size(r_hist,1);
a_hist = zeros(n_steps,1);
e_hist = zeros(n_steps,1);
i_hist = zeros(n_steps,1);
RAAN_hist = zeros(n_steps,1);
omega_hist = zeros(n_steps,1);
nu_hist = zeros(n_steps,1);
for k=1:n_steps
    [e_h, nu_h] = rv2coe(r_hist(k,:)', v_hist(k,:)', mu, e);
    a_hist(k) = a;
    e_hist(k) = e;
    i_hist(k) = i / deg;
    RAAN_hist(k) = RAAN / deg;
    omega_hist(k) =omega / deg;
    nu_hist(k) = nu_h / deg;
end

% MATLAB result output (plotting results) %
% 3D Trajectory with earth as refernce
figure('Name','3D Orbit','NumberTitle','off');
plot3(r_hist(:,1), r_hist(:,2), r_hist(:,3), '-b','LineWidth',1.6); hold on;
plot3(r_hist(1,1), r_hist(1,2), r_hist(1,3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
[XS,YS,ZS] = sphere(60);
surf(RE*XS, RE*YS, RE*ZS, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
text(r_hist(1,1), r_hist(1,2), r_hist(1,3),sprintf('  Orbital Period = %.5f min',T/60))
axis equal; grid on;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
title('3D Orbit Trajectory (ECI)');
% r(t) and v(t) history graphs
figure('Name','r and v histories','NumberTitle','off');
subplot(2,1,1);
plot(tt/60, r_mag, 'g','LineWidth',1.2);
xlabel('Time (min)'); ylabel('|r| (km)'); grid on;
title('Position Magnitude VS Time');
subplot(2,1,2);
plot(tt/60, v_mag, '-r','LineWidth',1.2);
xlabel('Time (min)'); ylabel('|v| (km/s)'); grid on;
title('Velocity Magnitude VS Time');
% Orbital elements time histories
figure('Name','Orbital Elements','NumberTitle','off');
subplot(3,2,1); plot(tt/60, a_hist,'-b'); xlabel('Time (min)'); ylabel('a (km)'); grid on;
subplot(3,2,2); plot(tt/60, e_hist,'-b'); xlabel('Time (min)'); ylabel('e'); grid on;
subplot(3,2,3); plot(tt/60, i_hist,'-b'); xlabel('Time (min)'); ylabel('i (deg)'); grid on;
subplot(3,2,4); plot(tt/60, RAAN_hist,'-b'); xlabel('Time (min)'); ylabel('\Omega (deg)'); grid on;
subplot(3,2,5); plot(tt/60, omega_hist,'-b'); xlabel('Time (min)'); ylabel('\omega (deg)'); grid on;
subplot(3,2,6); plot(tt/60, nu_hist,'-b'); xlabel('Time (min)'); ylabel('\theta (deg)'); grid on;
sgtitle('Time histories of classical orbital elements');

%% Transfer Orbit Propagator (Mission Phase 2)
% given data for mission phase
rp_parking = RE + 300;      % perigee radius = 300 km altitude [km]
i = 28.5 * deg;             % inclination [rad]
rp_GTO = RE + 300;          % Target perigee radius [km]
ra_GTO = 42164;             % Target apogee radius [km]

% use Hohmann transfer algorithem to find Δv at perigee for GTO
% 2. Semi-major axis of transfer ellipse
a_GTO = (rp_GTO + ra_GTO) / 2; %find Semi-major axis of transfer orbit (1)
v_circ = sqrt(mu / rp_parking); %Circular velocity in original orbit (2)
v_perigee_GTO = sqrt( mu * (2/rp_GTO - 1/a_GTO) ); %Velocity at perigee of the GTO (3)
dV_injection = v_perigee_GTO - v_circ; %Δv required at perigee (4)

% Algorithem two - propagate half revolution using two-body motion
% Define initial conditions (r, v) in perifocal plane
e_GTO = (ra_GTO - rp_GTO) / (ra_GTO + rp_GTO);
h = sqrt(mu * a_GTO * (1 - e_GTO^2));
r0 = [rp_GTO; 0; 0]; % start location = perigee
v0 = [0; h/rp_GTO; 0]; % start location = burn location
T = 2 * pi * sqrt(a_GTO^3 / mu); % Orbital period of transfer orbit
TOF = T / 2; %time of flight from perigee to apogee

% Propagate from perigee to apogee using ODE45 integration in 10e+4 steps
tspan = linspace(0, TOF, 100000);
y0 = [r0; v0];
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
[tt, yy] = ode45(@(t,y) twoBodyEOM(t,y,mu), tspan, y0, opts);

% get full velocity and position vector propagation history
r_hist = yy(:,1:3);
v_hist = yy(:,4:6);
r_mag = sqrt(sum(r_hist.^2,2));

% Plot Hohmann transfer trajectory in 3D
figure('Name','GTO Orbit Propagation','NumberTitle','off'); grid on;
plot3(r_hist(:,1), r_hist(:,2), r_hist(:,3), '-b'); hold on;
plot3(0,0,0,'ko','MarkerFaceColor','y','MarkerSize',8); % Earth center
plot3(r_hist(1,1), r_hist(1,2), r_hist(1,3),'go','MarkerFaceColor','g'); % Perigee
plot3(r_hist(end,1), r_hist(end,2), r_hist(end,3),'ro','MarkerFaceColor','r'); % Apogee
text(r_hist(1,1), r_hist(1,2),sprintf('  Δv required = %.5f Km/s ',dV_injection))
text(r_hist(end,1), r_hist(end,2),sprintf('  Orbital Period = %.5f hours \n   TOF = %.5f hours ',T/3600,TOF/3600))
xlabel('x (km)'); ylabel('y (km)');
title('GTO Trajectory (from Perigee to Apogee)');
legend('Trajectory','Earth','Perigee','Apogee');

%% Circularization and inclination reduction (Mission Phase 3)
% Given Phase data
a_GEO   = 42164.0;  % Semi-major axis of GEO (km)
e_GEO   = 0.0;      % Eccentricity of GEO
i_GEO   = 0.0;      % Target Inclination (deg)
r_perigee = RE + 300; % GTO Perigee Radius (km)
r_apogee  = a_GEO;      % GTO Apogee Radius (km) (Maneuver location)
i_GTO   = 28.5;     % Initial Inclination (deg, from previous problem)

% GTO and transfer calculated required data for manouver
a_GTO = (r_perigee + r_apogee) / 2; % GTO Semi-major axis
e_GTO = (r_apogee - r_perigee) / (r_apogee + r_perigee); % GTO Eccentricity
v_GTO_apo = sqrt(mu * (2/r_apogee - 1/a_GTO)); % GTO velocity at apogee
v_GEO = sqrt(mu / r_apogee); % velocity of target orbit

% Delta v calcilation for circulization and inclination reduction
delta_v_circularization = abs(v_GEO - v_GTO_apo); % delta
delta_i = deg2rad(i_GTO - i_GEO); 
delta_v_inclination = 2 * v_GEO * sin(delta_i / 2);
delta_v_combined = delta_v_inclination + delta_v_circularization;

% Plotting both manouver for visualization

% Define the GEO orbit trajectory 
theta_GEO = linspace(0, 2*pi, 100);
r_GEO_plot = a_GEO * [cos(theta_GEO); sin(theta_GEO); zeros(1, 100)];

% Rotate GEO plane to GTO plane for visualization clarity 
i_rad = deg2rad(i_GTO);
Omega_rad = deg2rad(270); % view angle
R_rot = [cos(Omega_rad) -sin(Omega_rad) 0; sin(Omega_rad) cos(Omega_rad) 0; 0 0 1] * ...
        [1 0 0; 0 cos(i_rad) sin(i_rad); 0 -sin(i_rad) cos(i_rad)];
r_GEO_plot_rot = R_rot * r_GEO_plot;

% Define the GTO Trajectory
E = linspace(-pi, pi, 100); % Eccentric Anomaly
r_GTO_magnitude = a_GTO * (1 - e_GTO * cos(E));
x_GTO_perifocal = a_GTO * (cos(E) - e_GTO);
y_GTO_perifocal = a_GTO * sqrt(1 - e_GTO^2) * sin(E);

% Find apogee and perigee location on trajectory
x_apo = a_GTO * (-1 - e_GTO); 
y_apo = 0;
x_peri = a_GTO * (1 - e_GTO);
y_peri = 0;
% generate 3d plot with found delta v
figure('Name', 'GTO to GEO Maneuver Visualization');
subtitle(sprintf('  Total Δv required for Phase 4 = %.5f Km/s ', delta_v_combined))
hold on;
plot3(x_GTO_perifocal, -zeros(size(y_GTO_perifocal)), y_GTO_perifocal, 'g--', 'LineWidth', 1, 'DisplayName', 'GTO Trajectory (28.5\circ Incl)');
[X_E, Y_E, Z_E] = sphere(60); % Plot Earth
surf(X_E * RE, Y_E * RE, Z_E * RE, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
plot3(r_GEO_plot(1, :), r_GEO_plot(2, :), r_GEO_plot(3, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Target GEO Orbit (0\circ Incl)'); % Plot the target GEO orbit (Equatorial)
plot3(-r_apogee, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Maneuver Point (Apogee)'); %plot Burn location
quiver3(-r_apogee, 0, 0, 0, 5000, 0, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'DisplayName', '\Delta V Burn'); % plot delta v vector
text(r_apogee, 0, sprintf('  Inclination reduction Δv required = %.5f Km/s ', delta_v_inclination))
text(x_apo, y_apo, sprintf('  Circularization Δv required = %.5f Km/s ', delta_v_circularization))
axis equal;
grid on;
title('GTO to GEO  Maneuver');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('Location', 'best');
view(15,15);
hold off;

%% functions used

function dydt = twoBodyEOM(~, y, mu)
% Two-body equations of motion for gravitational central force
r = y(1:3);
v = y(4:6);
rmag = sqrt(r'*r);
a = -mu/rmag^3 * r;
dydt = [v; a];
end

function [r_eci, v_eci] = coe2rv(a, e, i, RAAN, omega, nu, mu)
% algorithem one to Convert classical orbital elements to ECI position and velocity

p = a*(1 - e^2); % semi-latus rectum (1)
% position & velocity in perifocal frame (4)
r_pqw = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
v_pqw = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];
% Rotation matrices (5)
R3_W = [ cos(RAAN)  -sin(RAAN)   0;
    sin(RAAN)   cos(RAAN)   0;
    0           0       1];
R1_i = [1      0           0;
    0   cos(i)   -sin(i);
    0   sin(i)    cos(i)];
R3_w = [ cos(omega)   -sin(omega)    0;
    sin(omega)    cos(omega)    0;
    0              0       1];
% Transform PQW to ECI (6)
Q = R3_W * R1_i * R3_w;
r_eci = Q * r_pqw;
v_eci = Q * v_pqw;
end

function [ e_out, nu] = rv2coe(r, v, mu, e)
% Convert position & velocity vectors to get instantanous theta
% magnitudes
rmag = norm(r);
vmag = norm(v);
% eccentricity vector
evec = (1/mu) * ( (vmag^2 - mu/rmag)*r - (dot(r,v))*v );
e_out = e;
% true anomaly
if e > 1e-12
    nu = acos( dot(evec, r) / (e * rmag) );
    if dot(r, v) < 0
        nu = 2*pi - nu;
    end
end

end

