clc; clear; close all;

H = 21336;
[~, a, Pinf, rho_atm] = atmosisa(H, extended=true);

%% Fuel Properties
params.gamma = 1.3;       % ratio of specific heats
params.R_ = 287;          % gas constant
params.T0_ = 3440;        % flame temperature (?)
params.rho = 1855;        % propellant solid density
params.rho_g = 0.001;     % propellant gas density (zero or some small value)
params.a = 10.19*10^(-(3+6*0.433));  % burn rate termperature coefficient
params.n = 0.433;         % burn rate pressure exponent

%% Motor Geometry
params.A_t = 0.00493;     % nozzle throat area
params.r_max = 0.127;     % Maximum radius
params.V_hardware = 0.2*pi*params.r_max^2; % extra hardware volume, if any

%% Grain Geometry
params.len = 1;          % motor grain length, etc.
params.cs = finocyl_pointy(4, 5, 4, 0.4);

%% something
P0      = 100;           % [Pa], for example

X0      = [ P0; ]; % combine into the state vector
tspan   = [0, 30];        % want to simulate 2 seconds

[tSol, XSol] = ode45(@(t,X) RocketSystem_Tapered(t, X, params), tspan, X0, odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Events', @(t,X) EVENT_ACCELPHASECOMPLETE(t, X, params.r_max)));

P_sol    = XSol(:,1);   % pressure over time
r_in_sol = XSol(:,2);   % inner radius over time
r_out_sol = XSol(:,3);  % outer radius over time
len_sol = XSol(:,4);    % length over time


P0      = P_sol(end);           % [Pa], for example
len0    = 0.2;
X0      = [ P0; len0 ]; % combine into the state vector
tspan   = [tSol(end), tSol(end)+4];        % want to simulate 2 seconds

[tSol2, XSol2] = ode45(@(t,X) RocketSystem_Cylinder(t, X, params), tspan, X0, odeset('Events', @(t,X) EVENT_CRUISEPHASECOMPLETE(t, X, len0)));

P_sol   = [P_sol; XSol2(:,1)];
tSol_tot   = [tSol; tSol2];

Ab_sol = 2 * pi * r_in_sol .* len_sol;
t_angle = atan(len_sol./(r_out_sol - r_in_sol)).*(180/pi);          % taper angle

figure;
plot(tSol_tot, P_sol, 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 12);
ylabel('Pressure [Pa]', 'FontSize', 12);
title('Chamber Pressure vs Time', 'FontSize', 12);
grid on;

figure;
% --- Top subplot: Pressure ---
subplot(4,1,1);
plot(tSol, r_in_sol, 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 12);
ylabel('r_{in} [m]', 'FontSize', 12);
title('Grain Inner Radius vs Time', 'FontSize', 12);
grid on;

% --- Bottom subplot: Inner radius ---
subplot(4,1,2);
plot(tSol, r_out_sol, 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 12);
ylabel('r_{in} [m]', 'FontSize', 12);
title('Grain Outer Radius vs Time', 'FontSize', 12);
grid on;

subplot(4,1,3);
plot(tSol, len_sol, 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 12);
ylabel('r_{in} [m]', 'FontSize', 12);
title('Grain Length vs Time', 'FontSize', 12);
grid on;

subplot(4,1,4);
plot(tSol(1:end-1), t_angle(1:end-1), 'LineWidth', 2);
xlabel('Time [s]', 'FontSize', 12);
ylabel('r_{in} [m]', 'FontSize', 12);
title('Taper Angle vs Time', 'FontSize', 12);
grid on;

N = length(r_in_sol);

