clear all;
close all;
clc;

%% =========================================================
%  Bonus Point: CR3BP Transfer
% ==========================================================

%% =========================================================
%  PARAMETERS
% ==========================================================
mu_M = 4902.8;        % Moon gravitational parameter [km^3/s^2]
R_M  = 1737.4;        % Moon mean radius [km]
h1   = 100;           % LM orbit altitude [km]
h2   = 400;           % MS orbit altitude [km]
r_LM = R_M + h1;      % LM orbital radius [km]
r_MS = R_M + h2;      % MS orbital radius [km]
phi0 = 10 * (pi/180); % Initial phasing angle [rad]
mu = 0.01215;         % Mass Ratio Earth-Moon system
d_EM = 384400;        % Earth–Moon distance [km]

%% =========================================================
%  POINT 1 – Hohmann Transfer: Delta-V and Time of Flight
% ==========================================================
v_LM = sqrt(mu_M / r_LM);   % Circular velocity on LM orbit [km/s]
v_MS = sqrt(mu_M / r_MS);   % Circular velocity on MS orbit [km/s]

a = (r_LM + r_MS) / 2;                   % Transfer ellipse SMA [km]
v_perigee = sqrt(mu_M * (2/r_LM - 1/a)); % Velocity at perigee [km/s]
v_apogee  = sqrt(mu_M * (2/r_MS - 1/a)); % Velocity at apogee  [km/s]

delta_V1      = v_perigee - v_LM;        % First impulse  [km/s]
delta_V2      = v_MS - v_apogee;         % Second impulse [km/s]
total_delta_V = delta_V1 + delta_V2;     % Total delta-V  [km/s]
delta_t = pi * sqrt(a^3 / mu_M);         % Transfer time of flight [s]

fprintf('===== POINT 1 =====\n');
fprintf('Delta_V1          = %.6f km/s\n', delta_V1);
fprintf('Delta_V2          = %.6f km/s\n', delta_V2);
fprintf('Total Delta-V     = %.6f km/s\n', total_delta_V);
fprintf('Time of flight    = %.4f s  (%.2f h)\n\n', delta_t, delta_t/3600);

%% =========================================================
%  POINT 2 – Required Phase Angle and Waiting Time
% ==========================================================
n_LM = sqrt(mu_M / r_LM^3);  % Mean motion LM [rad/s]
n_MS = sqrt(mu_M / r_MS^3);  % Mean motion MS [rad/s]

% Angle swept by MS during transfer
theta_MS_transfer = n_MS * delta_t;
% For rendezvous
phi_required = pi - theta_MS_transfer; % [rad]

% Gap to close by waiting on the initial orbit
delta_phi = mod(phi0 - phi_required, 2*pi); % [rad]
n_rel  = n_LM - n_MS;                       % Relative angular velocity [rad/s]
t_wait = delta_phi / n_rel;                 % Waiting time [s]

fprintf('===== POINT 2 =====\n');
fprintf('Required phase angle = %.6f deg\n', phi_required * 180/pi);
fprintf('Initial phase angle  = %.6f deg\n', phi0 * 180/pi);
fprintf('Delta-phi to cover   = %.6f deg\n', delta_phi * 180/pi);
fprintf('Waiting time         = %.4f s  (%.2f h)\n\n', t_wait, t_wait/3600);


%% =========================================================
%  POINT 3 – Orbit Plot
% ==========================================================

figure('Name','Point 3 – Hohmann Transfer','Color','w');
hold on; grid on; axis equal;

theta_full = linspace(0, 2*pi, 400);
R_M7=R_M*0.7;
fill(R_M7*cos(theta_full), R_M7*sin(theta_full), [0.75 0.75 0.75], ...
    'EdgeColor','k','DisplayName','Moon');

plot(r_LM*cos(theta_full), r_LM*sin(theta_full), 'b--', ...
    'LineWidth',0.9,'DisplayName','LM Orbit  (h_1=100 km)');
plot(r_MS*cos(theta_full), r_MS*sin(theta_full), 'r--', ...
    'LineWidth',0.9,'DisplayName','MS Orbit  (h_2=400 km)');

ecc = (r_MS - r_LM)/(r_MS + r_LM);
p   = a*(1 - ecc^2);
th_tr = linspace(0, pi, 300);
r_tr  = p ./ (1 + ecc*cos(th_tr));
plot(r_tr.*cos(th_tr), r_tr.*sin(th_tr), 'g-', ...
    'LineWidth',1.8,'DisplayName','Hohmann Transfer Ellipse');
plot(r_LM, 0, 'bo','MarkerFaceColor','b','MarkerSize',8, ...
    'DisplayName','LM – Transfer Departure');
plot(r_MS*cos(phi_required), r_MS*sin(phi_required), 'ro', ...
    'MarkerFaceColor','r','MarkerSize',8,'DisplayName','MS – Transfer Departure');
plot(-r_MS, 0, 'b*','MarkerSize',10,'LineWidth',1.5, ...
    'DisplayName','Rendezvous Location');

% phi_required arc
r_arc      = 0.68 * r_MS;
theta_arc  = linspace(0, phi_required, 200);
plot([0, r_arc], [0, 0], 'm:', 'LineWidth', 1.2, 'HandleVisibility','off');
plot([0, r_arc*cos(phi_required)], [0, r_arc*sin(phi_required)], 'm:', ...
    'LineWidth', 1.2, 'HandleVisibility','off');
plot(r_arc*cos(theta_arc), r_arc*sin(theta_arc), 'm-', ...
    'LineWidth', 2.0, 'DisplayName', ...
    sprintf('\\phi_{req} = %.1f°', phi_required*180/pi));

% Arrowhead and Label
arrow_len = 0.035 * r_MS;
tang      = phi_required + pi/2;  
tip       = [r_arc*cos(phi_required), r_arc*sin(phi_required)];
p1        = tip - arrow_len * [cos(tang - 0.4), sin(tang - 0.4)];
p2        = tip - arrow_len * [cos(tang + 0.4), sin(tang + 0.4)];
fill([tip(1), p1(1), p2(1)], [tip(2), p1(2), p2(2)], 'm', ...
    'EdgeColor','m', 'HandleVisibility','off');
theta_mid = phi_required / 2;
r_label   = r_arc + 130;
text(r_label*cos(theta_mid), r_label*sin(theta_mid), ...
    sprintf('\\phi_{req} = %.1f°', phi_required*180/pi), ...
    'Color','m', 'FontSize', 9, 'FontWeight','bold', ...
    'HorizontalAlignment','center');

xlabel('X [km]','FontWeight','bold'); ylabel('Y [km]','FontWeight','bold');
title('Hohmann Transfer Geometry - Three Body Problem','FontSize',13);
legend('Location','northeastoutside','FontSize',8);
xlim([-(r_MS+500),(r_MS+500)]); ylim([-(r_MS+500),(r_MS+500)]);
hold off;

%% =========================================================
%  POINT 3 – SHOOTING METHOD IN THE CR3BP
% ==========================================================

% Non-dimensionalization
MU_tot = mu_M / mu;       % Total gravitation parameter of the system [km^3/s^2]
DU = d_EM;                % Distance Unit
TU = sqrt(DU^3 / MU_tot); % Time Unit
VU = DU/TU;               % Velocity Unit

r_LM_nd = r_LM/DU;        % Dimensionless radius of parking orbit (LL0)
r_MS_nd = r_MS/DU;        % Dimensionless radius of the MS orbit
dt_guess = delta_t/TU;    % Dimensionless Hohmann time of flight
v_per_nd = v_perigee/VU;  % Dimensionless velocity at the perigee of the Hohmann ellipse

% Rotating Reference Frame (Synodic System)
% We set the intial state of the S/C at deperture
x0_nd = (1-mu)+r_LM_nd; % X Position: Moon's center + LLO radius
y0_nd = 0;              % Y Position: starting exactly on the X-radius
z0_nd = 0;              % Z Position: the motion is 2D, so out-of-plane altitude is zero

vxo_nd = 0;                      % X Velocity: at the perigee, the radial velocity is zero
vy0_guess_nd = v_per_nd-r_LM_nd; % Y Velocity: converted to the rotating reference system
vz0_nd = 0;                      % Z Velocity: zero for planar motion

% The provisional initial state vector
X0_guess = [x0_nd; y0_nd; z0_nd; vxo_nd; vy0_guess_nd; vz0_nd];

% Define Target
x_target_nd = (1-mu) - r_MS_nd; % where the LM needs to go

% Shooting Algorithm (fsolve Solver)
% The two variables that fsolve needs to 'guess' and correct
guess_vars = [vy0_guess_nd; dt_guess]; %[Initial Y Velocity; Time of Flight]

% Set strict mathematical tolerances (1e-10) for the target error
options_fsolve = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-10, 'TolX', 1e-10);

% You can't pass the shooting_method function directly to fsolve because 
% it wouldn't know which parameteres to vary, so by doing so, it tells it
% to vary only 'vars' each time it calls the shooting_method
target_fun = @(vars) shooting_method(vars, X0_guess, x_target_nd, mu);

% Run the solver: fsolve will iterate untile the error goes to zero
[sol_vars, fval, exitflag] = fsolve(target_fun, guess_vars, options_fsolve);

% Exctract the winning solution from the sol_vars vector
vy0_cr3bp = sol_vars(1);   % The perfect Y velocity to compensate for Earth's gravity
dt_cr3bp_nd = sol_vars(2); % The real time of flight in the CR3BP

% Prepare the final state vector by overwritign the guessed
X0_cr3bp = X0_guess;     
X0_cr3bp(5) = vy0_cr3bp; % Index 5 corresponds to vy

% Integrate the EoM with high tolerences
options_ode = odeset('RelTol',1e-11, 'AbsTol',1e-11);
[t_cr3bp, X_cr3bp_nd] = ode45(@(t,X) cr3bp(t,X,mu),[0,dt_cr3bp_nd], X0_cr3bp, options_ode);

% --- Delta V Calculation (Maneuver Costs) ---
% Maneuver 1: Departure from LLO
v1_rot_nd = (v_LM/VU)-r_LM_nd; % Circular LLO velocity in the rotating frame
dV1_cr3bp = abs(vy0_cr3bp - v1_rot_nd)*VU; % Delta V: velocity difference [Km/s]

% Maneuver 2: Arrival at the MS orbit
v2_rot_nd = -(v_MS/VU)+r_MS_nd; % Circular MS velocity in the rotating frame
vy_arrival_cr3bp = X_cr3bp_nd(end,5); % Actual velocity upon hitting the target 
dV2_cr3bp = abs(vy_arrival_cr3bp-v2_rot_nd)*VU; % Braking Delta V [Km/s]

% Integrate the EoM without shooting method
[t_cr3bp, X_inco_nd] = ode45(@(t,X) cr3bp(t,X,mu),[0,dt_guess], X0_guess, options_ode);

% Incorrect final trajectory position
x_end_nd = X_inco_nd(end,1);
y_end_nd = X_inco_nd(end,2);

% Calculate the distance from the center of the Moon
r_final_nd = sqrt((x_end_nd-(1-mu))^2+y_end_nd^2); 
r_final_km = r_final_nd*DU;  %[km]

% Radial error with respect to the target orbit
error_radius_km = abs(r_final_km - r_MS);

% Print final results to the command window
fprintf('Time of flight CR3BP: %.2f s (vs %.2f s Hohmann 2-Body)\n', dt_cr3bp_nd*TU, delta_t);
fprintf('Delta_V1 CR3BP: %.4f km/s (vs %.4f km/s Hohmann)\n', dV1_cr3bp, delta_V1);
fprintf('Delta_V2 CR3BP: %.4f km/s (vs %.4f km/s Hohmann)\n', dV2_cr3bp, delta_V2);
fprintf('Total Delta_V: %.4f km/s (vs %.4f km/s Hohmann)\n', dV1_cr3bp+dV2_cr3bp, total_delta_V);
fprintf('Target Orbit Radius:     %.3f km\n', r_MS);
fprintf('Achieved Radius (CR3BP):  %.3f km\n', r_final_km);
fprintf('Radial Error: %.4f m\n\n', error_radius_km*1000);

% --- Trajectory Plotting ---
figure(1); hold on;

% Shift the origin from the Earth-Moon barycenter to the Moon's center and
% re-multiply bu DU
x_plot_cr3bp = (X_cr3bp_nd(:,1) - (1 - mu)) * DU;
y_plot_cr3bp = X_cr3bp_nd(:,2) * DU;
plot(x_plot_cr3bp, y_plot_cr3bp, 'y--', 'LineWidth', 2, 'DisplayName', 'CR3BP Real Trajectory');
legend('Location', 'best');

%% =========================================================
%  POINT 4 – Calculating the Jacobi Constant
% ==========================================================

% Exctract component in rotating reference frame
x_rot = X_cr3bp_nd(:,1);
y_rot = X_cr3bp_nd(:,2);
z_rot = X_cr3bp_nd(:,3);
vx_rot = X_cr3bp_nd(:,4);
vy_rot = X_cr3bp_nd(:,5);
vz_rot = X_cr3bp_nd(:,6);

N = length(t_cr3bp); % Total number of integrations steps calculated by ode45
C = zeros(N,1);      % Pre-allocate memory for the Jacobi Constant array

for i = 1:N

   r1_rot = sqrt((mu+x_rot(i))^2 + y_rot(i)^2 + z_rot(i)^2);     % Distance from the first primary (Earth)
   r2_rot = sqrt((x_rot(i)-(1-mu))^2 + y_rot(i)^2 + z_rot(i)^2); % Distance from the second primary (Moon)

   V = vx_rot(i)^2+vy_rot(i)^2+vz_rot(i)^2; % Squared velocity magnitude in the rotating reference frame (V^2)
   U = 0.5*(x_rot(i)^2+y_rot(i)^2) + (1-mu)/r1_rot + mu/r2_rot; % Effective pseudo-potentional (U)
   C(i) = 2*U - V; % Jacobi Constant 
end

% Calculate the absolute error relative to the initial energy state (t=0)
C_error = abs(C-C(1));

% Plotting the Jacobi Constant Error
figure('Name', 'Jacobi Error Constant', 'Color', 'w', 'Position', [100, 100, 900, 700]);
plot(t_cr3bp*TU, C_error,'r-','LineWidth',1.5);
title('Jacobi constant convertion error ');
xlabel('Time [s]');
ylabel('Absolute Error |C(t)-C(0)|');
grid on;

% Extract and print the maxium error value to validate the solver's tolerances
max_err = max(C_error);
fprintf('Max Error of Jacobi Constant: %.2e', max_err);

%% =========================================================
%  LOCAL FUNCTIONS
% ==========================================================

% CR3BP Equations of Motion in the rotating (synodic) reference frame
% This function calculates the state derivative vector for ode45
function dX = cr3bp(~,X,mu)

% Extract current position and velocity from the state vector
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);

% --- Distances to the Primaries ---
r1 = sqrt((mu+x)^2 + y^2 + z^2); %r1: Distance from the spacecraft to the first primary (Earth, at -mu)

r2 = sqrt((x-(1-mu))^2 + y^2 + z^2); % r2: Distance from the spacecraft to the second primary (Moon, at 1-mu)

% --- Partial Derivatives of the Effective Potential (U) ---
% These terms include both the gravitational attraction of the two masses 
% and the centrifugal force caused by the rotating frame.
dU_dx = x - (1 - mu)*(x + mu)/r1^3 - mu*(x - (1 - mu))/r2^3;
dU_dy = y - (1 - mu)*y/r1^3 - mu*y/r2^3;
dU_dz = -(1 - mu)*z/r1^3 - mu*z/r2^3;

% --- Accelerations ---
% The terms 2*vy and -2*vx represent the Coriolis acceleration, 
% which acts purely in the rotating frame and is velocity-dependent.
ax = dU_dx + 2*vy;
ay = dU_dy - 2*vx;
az = dU_dz;

% Assemble and return the state derivative vector [Velocity; Acceleration]
dX = [vx; vy; vz; ax; ay; az];

end


% Target function for the fsolve algorithm.
% It propagates a test trajectory and returns the miss distance (error).
function err = shooting_method(vars, X0, x_target, mu)

% 1. Extract the trial variables generated by the solver
vy0 = vars(1); % Trial initial velocity along the y-axis
tof = vars(2); % Trial time of flight

% 2. Update the initial state vector with the trial velocity 
X0(5) = vy0;

% 3. The "Shoot" (Numerical Integration)
% Propagate the spacecraft using the trial velocity and time
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
[~, X_out] = ode45(@(t,X) cr3bp(t,X,mu), [0, tof], X0, options);

% 4. Extract the final position at the end of the time of fligh 
x_final = X_out(end, 1);
y_final = X_out(end, 2);
    
% 5. Evaluate the Error
% We want the spacecraft to cross the x-axis (y = 0) exactly at x_target.
% fsolve will continuously adjust 'vars' until this vector equals [0; 0].
err = [x_final - x_target; 
           y_final];
end
