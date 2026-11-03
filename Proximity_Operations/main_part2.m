clear all;
close all;
clc;

%% =========================================================
%  PART 2 – Relative Motion and HCW Dynamics
%
%  The LM is assumed to end the Hohmann transfer with residual
%  relative position and velocity errors with respect to the MS.
%
%  The relative dynamics are modeled using the linearized
%  Hill–Clohessy–Wiltshire (HCW) equations in the LVLH frame:
%
%     x-axis -> along-track direction (V-bar)
%     y-axis -> out-of-plane direction (H-bar)
%     z-axis -> radial direction toward the Moon (R-bar)
% ==========================================================

% ==========================================================
%   POINT 1 – Random Initial Errors
% ==========================================================

% --- Relative position error (random direction, constrained norm) ---
delta_r_dir = randn(3,1);
delta_r_dir = delta_r_dir/norm(delta_r_dir);   % Normalization 

delta_r_min = 0.5;       % 500 m = 0.5 km
delta_r_max = 1;         % 1000 m = 1 km
rho_mag = delta_r_min + (delta_r_max-delta_r_min)*rand;  % Random position error magnitude inside required range [km]

delta_r0 = rho_mag * delta_r_dir;   % Initial relative position error vector [km]

% --- Relative velocity error (random direction, constrained norm) ---
delta_v_dir = randn(3,1);
delta_v_dir = delta_v_dir/norm(delta_v_dir);    % Normalization 

delta_v_min = 0.0001;    % 0.1 m/s = 1e-4 km/s
delta_v_max = 0.001;     % 1   m/s = 1e-3 km/s
v_mag = delta_v_min + (delta_v_max-delta_v_min)*rand;    % Random velocity error magnitude inside required range [km/s]

delta_v0 = v_mag * delta_v_dir;     % Initial relative velocity error vector [km/s]

% --- Initial relative state vector (LVLH frame) ---
% State = [x; y; z; x_dot; y_dot; z_dot]
s0 = [delta_r0; delta_v0];

% Reference orbit = MS circular lunar orbit
mu_M = 4902.8;            % Moon gravitational parameter [km^3/s^2]
R_M  = 1737.4;            % Moon mean radius [km]
h2   = 400;               % MS orbit altitude [km]
r_MS = R_M + h2;          % MS orbital radius [km]
n = sqrt(mu_M / r_MS^3);  % Mean motion MS [rad/s]


% =========================================================
%  HCW Analytical Propagation using STM
%  Relative state evolution:  s(t) = Phi(t) * s0
%  where Phi(t) is the HCW State Transition Matrix.
% =========================================================

% initial relative position -> relative position
phi_rr = @(t) [1, 0, 6*n*t-6*sin(n*t);
    0, cos(n*t), 0;
    0, 0, 4-3*cos(n*t)];

% initial relative velocity -> relative position
phi_rv = @(t) [((4*sin(n*t))/n)-3*t, 0, (2*(1-cos(n*t)))/n;
    0, sin(n*t)/n, 0;
    (2*(cos(n*t)-1))/n, 0, sin(n*t)/n];

% initial relative position -> relative velocity
phi_vr = @(t) [0, 0, 6*n*(1-cos(n*t));
    0, -n*sin(n*t), 0;
    0, 0, 3*n*sin(n*t)];

% initial relative velocity -> relative velocity
phi_vv = @(t) [4*cos(n*t)-3, 0, 2*sin(n*t);
    0, cos(n*t), 0;
    -2*sin(n*t), 0, cos(n*t)];

% ==========================================================
%   POINT 2 – Propagate the relative motion
% ==========================================================

T = 2*pi / n;                 % MS orbital period [s]
t_vec = linspace(0,3*T,2000); % Time span vector for propagation
s = zeros(length(t_vec),6);   % State matrix inizialization  

for i = 1:length(t_vec)
    t = t_vec(i);

% Complete HCW State Transition Matrix (6x6)
phi = [phi_rr(t), phi_rv(t);
    phi_vr(t), phi_vv(t)];

s(i,:) = phi*s0;              % Propagated relative state at current time
end

% =========================================================
%  Free-Drift Relative Motion Plot
% ==========================================================

figure('Name','PART 2 – HCW Free Drift','Color','w');

hold on;
grid on;
axis equal;

% 3D free-drift trajectory
plot3(s(:,1)*1000, s(:,2)*1000, s(:,3)*1000, 'b-', 'LineWidth',1.8, 'DisplayName','LM Relative Trajectory');
% Initial position
plot3(s(1,1)*1000, s(1,2)*1000, s(1,3)*1000, 'go', 'MarkerFaceColor','g', 'MarkerSize',8, 'DisplayName','Initial Position');
% Final position
plot3(s(end,1)*1000, s(end,2)*1000, s(end,3)*1000, 'ro', 'MarkerFaceColor','r', 'MarkerSize',8, 'DisplayName','Final Position');

 
xlabel('x - Along-track (v-bar) [m]','FontWeight','bold');
ylabel('y - Out-of-plane [m]','FontWeight','bold');
zlabel('z - Radial [m]','FontWeight','bold');

title('HCW Free-Drift Relative Motion','FontSize',13);

legend('Location','best');

view(35,25);

hold off;

%% =========================================================
%   POINT 3 – Two-impulse corrective maneuver
% ==========================================================

% =========================================================
% Design of a two-impulse maneuver using the HCW framework
% to target a 50 m V-bar hold point ahead of the MS.
% =========================================================

% Set the Target Conditions
r_final = [0.05,0,0]'; % Position vector in LVLH frame [km]
v_final = [0,0,0]';    % Final velocity vector for hold position [km/s]

% Transfer Time
t_man = T/4;           % A quarter of the orbital period

% =========================================================
% Compute the First Impulse (Departure Burn)
% From: r_final = phi_rr*r_0 + phi_rv*(v_0+delta_v1)
% This impulse puts the LM on a collision course with the 
% target point.
% =========================================================

% First Impulse
delta_v1 = phi_rv(t_man) \ (r_final-phi_rr(t_man)*delta_r0)-delta_v0;
delta_v1_mag = norm(delta_v1);

% =========================================================
% Compute the Second Impulse (Arrival Burn)
% This impulse required to stop the LM 
% =========================================================

% Second Impulse 
delta_v2 = -(phi_vr(t_man)*delta_r0 + phi_vv(t_man)*(delta_v0 + delta_v1));
delta_v2_mag = norm(delta_v2);

% Total Cost
total_delta_v = norm(delta_v1) + norm(delta_v2); % Sum of magnitude of the two burns.

fprintf('===== POINT 3 =====\n');
fprintf('Delta_V1          = %.6f m/s\n', delta_v1_mag*1000);
fprintf('Delta_V2          = %.6f m/s\n', delta_v2_mag*1000);
fprintf('Total Delta-V     = %.6f m/s\n\n', total_delta_v*1000);


%% =========================================================
%  POINT 4 – Straight-Line V-bar Approach
% =========================================================
%
%  In this phase, the LM performs the final approach from the hold point
%  (computed in Point 3) to the docking port of the MS, located at the origin
%  of the LVLH frame.
%
%  The maneuver is designed as a forced straight-line motion along the V-bar
%  (along-track direction, x-axis), using a sequence of impulsive corrections.
%
%  At each time step:
%    - a desired velocity is computed to ensure linear convergence to the target
%    - an impulsive correction ΔV is applied
%    - the relative motion is propagated using HCW equations over Δt_imp
%
%  The total ΔV is obtained by summing all impulsive corrections.
% ==========================================================


% Initial condition = hold point from Point 3
r_curr = r_final;      % [km]
v_curr = v_final;      % should be 0

% Target docking point
r_target = [0;0;0];

% Control parameters
N_imp = 20;         % number of impulses
dt_imp = 30;        % time step [s]

traj_vbar = zeros(N_imp,3);

DeltaV_vbar = 0;

for k = 1:N_imp

    t_remain = (N_imp - k + 1)*dt_imp;      % Remaining time
    v_des = [(r_target(1)-r_curr(1))/t_remain; 0; 0];    % Desired velocity
    dv = v_des - v_curr;                    % Required impulse

    DeltaV_vbar = DeltaV_vbar + norm(dv);   % Total DeltaV computation

    v_curr = v_curr + dv;                   % Impulse apply

    % --- HCW propagation for dt_imp ---

    phi_rr = [1, 0, 6*n*dt_imp-6*sin(n*dt_imp);
        0, cos(n*dt_imp), 0;
        0, 0, 4-3*cos(n*dt_imp)];

    phi_rv = [((4*sin(n*dt_imp))/n)-3*dt_imp, 0, (2*(1-cos(n*dt_imp)))/n;
        0, sin(n*dt_imp)/n, 0;
        (2*(cos(n*dt_imp)-1))/n, 0, sin(n*dt_imp)/n];

    phi_vr = [0, 0, 6*n*(1-cos(n*dt_imp));
        0, -n*sin(n*dt_imp), 0;
        0, 0, 3*n*sin(n*dt_imp)];

    phi_vv = [4*cos(n*dt_imp)-3, 0, 2*sin(n*dt_imp);
        0, cos(n*dt_imp), 0;
        -2*sin(n*dt_imp), 0, cos(n*dt_imp)];

    % Propagation
    r_next = phi_rr*r_curr + phi_rv*v_curr;
    v_next = phi_vr*r_curr + phi_vv*v_curr;

    traj_vbar(k,:) = r_next';   % Trajectory State for plotting

    % Update state
    r_curr = r_next;
    v_curr = v_next;
end

fprintf('===== POINT 4 =====\n');
fprintf('Total Delta-V (V-bar approach) = %.6f m/s\n\n', DeltaV_vbar*1000);

% ======= TRAJECTORY PLOT =====================================

figure('Name','V-bar Approach 3D','Color','w');
hold on; grid on; axis equal;

plot3(traj_vbar(:,1)*1000, traj_vbar(:,2)*1000, traj_vbar(:,3)*1000, ...
    'b-o','LineWidth',1.5,'DisplayName','Trajectory');

% Docking port (origine)
plot3(0,0,0,'rs','MarkerFaceColor','r','MarkerSize',8,...
    'DisplayName','Docking Port');

% Starting Point (hold point)
plot3(traj_vbar(1,1)*1000, traj_vbar(1,2)*1000, traj_vbar(1,3)*1000, ...
    'go','MarkerFaceColor','g','MarkerSize',8,...
    'DisplayName','Hold Point');

xlabel('x (along-track) [m]','FontWeight','bold');
ylabel('y (out-of-plane) [m]','FontWeight','bold');
zlabel('z (radial) [m]','FontWeight','bold');

title('3D Straight-Line V-bar Approach');

legend('Location','best');

view(40,25);


%% ========================================================
%  POINT 5.2 – Full Nonlinear Relative Motion (LVLH Frame)
% =========================================================
%
%  We propagate the relative motion using the exact nonlinear equations
%  derived from the transformation between inertial and rotating frames.
%
%  Starting from:
%
%     d(rho)/dt_inertial = drho + omega x rho
%
%     d^2(rho)/dt^2_inertial =
%           d^2rho/dt^2
%         + 2 omega x drho/dt
%         + omega x (omega x rho)
%
%  Substituting the gravitational acceleration:
%
%     d^2(rho)/dt^2_inertial =
%        -mu * ( r_chaser/|r_chaser|^3 - r_target/|r_target|^3 )
%
%  we obtain the FULL NONLINEAR RELATIVE MOTION equation:
%
%     rho_ddot =
%        - mu * ( r_chaser/|r_chaser|^3 ) 
%        + mu * ( r_target/|r_target|^3 ) 
%        - 2 omega x rho_dot 
%        - omega x (omega x rho)
%
%  This model does NOT rely on linearization (unlike HCW).
% ==========================================================

omega = [0; -n; 0];              % LVLH angular velocity [rad/s]
r_target  = [0; 0; -r_MS];       % MS position (z = radial axis) [km]

f_nonlin = @(~,s) [ s(4:6);
    - mu_M*( r_target + s(1:3) ) / norm(r_target + s(1:3))^3 ...
    + mu_M*r_target/norm(r_target)^3 ...
    - 2*cross(omega, s(4:6)) ...
    - cross(omega, cross(omega, s(1:3))) ...
    ];

opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[t_nonlin, s_nonlin] = ode45(f_nonlin, t_vec, s0, opts);
r_nonlin = s_nonlin(:,1:3);

% =========================================================
%  ERROR NON LINEAR - HCW
% =========================================================

r_hcw = s(:,1:3);                   % HCW relative position [km]
err = vecnorm(r_nonlin - r_hcw,2,2);    % position error at each time [km]

fprintf('===== POINT 5.2 =====\n');
fprintf('Max Error Non Linear - HCW = %.6f m\n\n', max(err)*1000);

% ======= TRAJECTORY PLOT =================================

figure('Name','HCW - Full Nonlinear Relative Motion','Color','w');
hold on; grid on; axis equal;

plot3(r_nonlin(:,1)*1000, r_nonlin(:,2)*1000, r_nonlin(:,3)*1000, ...
    'r-','LineWidth',1.5,'DisplayName','Full nonlinear (general)');

plot3(r_hcw(:,1)*1000, r_hcw(:,2)*1000, r_hcw(:,3)*1000, ...
    'b--','LineWidth',1.5,'DisplayName','HCW');

% Initial position
plot3(s(1,1)*1000, s(1,2)*1000, s(1,3)*1000, 'go', 'MarkerFaceColor','g', 'MarkerSize',8, 'DisplayName','Initial Position');
% Final position
plot3(s(end,1)*1000, s(end,2)*1000, s(end,3)*1000, 'ro', 'MarkerFaceColor','r', 'MarkerSize',8, 'DisplayName','Final Position');

xlabel('x (along-track) [m]');
ylabel('y (out-of-plane) [m]');
zlabel('z (radial) [m]');

title('HCW - Full Nonlinear Relative Motion');
legend('Location','best');

view(40,25);

% ======= ERROR PLOT ======================================

figure('Name','Error Non Linear - HCW','Color','w');
plot(t_vec/3600, err*1000,'k','LineWidth',1.5);
grid on;

xlabel('Time [h]');
ylabel('Error [m]');
title('HCW Deviation from Full Nonlinear Model');

%% =========================================================
%   POINT 5.3 – Nonlinear Two-Impulse Maneuver
% =========================================================
%
%  The first impulse ΔV1 is computed iteratively.
%  An initial guess is corrected using fsolve, while the trajectory
%  is propagated with ode45, until the final position matches the target.
%
%  The second impulse ΔV2 is then applied to cancel the residual velocity.
% =========================================================

r_target = r_final;   % [km]
t_man = T/4;          % [s]

% Function: position error at final time
fun = @(dv1) final_position_error(dv1, delta_r0, delta_v0, t_man, f_nonlin, opts, r_target);
dv1_guess = delta_v1;   % Initial guess (using HCW solution as starting point)

% Solve nonlinear targeting problem
dv1_NL = fsolve(fun, dv1_guess);

% Propagate with optimal dv1_NL
s_init = [delta_r0; delta_v0 + dv1_NL];
[~, s_temp] = ode45(f_nonlin, [0 t_man], s_init, opts);
v_final_NL = s_temp(end,4:6)';

dv2_NL = -v_final_NL;   % Second impulse (to stop)

% =========================================================
%  COST AND COMPARISON
% ==========================================================

% Total Delta-V (nonlinear)
DV_tot_NL = norm(dv1_NL) + norm(dv2_NL);

% HCW reference
DV_tot_HCW = norm(delta_v1) + norm(delta_v2);

% Errors
err_dv1 = norm(dv1_NL - delta_v1);
err_dv2 = norm(dv2_NL - delta_v2);
err_DV_tot = abs(DV_tot_NL - DV_tot_HCW);

fprintf('===== POINT 5.3 (NONLINEAR) =====\n');
fprintf('Total Delta-V (NL)   = %.6f m/s\n', DV_tot_NL*1000);
fprintf('Total Delta-V (HCW)  = %.6f m/s\n', DV_tot_HCW*1000);
fprintf('DeltaV total error   = %.6f m/s\n\n', err_DV_tot*1000);

fprintf('DeltaV1 difference   = %.6f m/s\n', err_dv1*1000);
fprintf('DeltaV2 difference   = %.6f m/s\n\n', err_dv2*1000);

% ======= ERROR BAR PLOT =====================================

figure('Name','DeltaV Error HCW vs Nonlinear','Color','w');

bar([err_dv1, err_dv2, err_DV_tot]*1000);
grid on;

xticklabels({'DeltaV1','DeltaV2','Total DeltaV'});
ylabel('Difference [m/s]');
title('HCW vs Nonlinear Impulse Difference');

%% =========================================================
%   POINT 5.4 – Nonlinear V-bar approach
% ==========================================================

% --- Initial condition (hold point) ---
r_curr_NL = r_final;
v_curr_NL = v_final;

r_target = [0; 0; 0];

traj_vbar_NL = zeros(N_imp,3);
DeltaV_vbar_NL = 0;

for k = 1:N_imp

    t_remain = (N_imp - k + 1)*dt_imp;

    % Desired velocity (same logic as HCW)
    v_des = [(r_target(1)-r_curr_NL(1))/t_remain; 0; 0];

    dv = v_des - v_curr_NL;

    DeltaV_vbar_NL = DeltaV_vbar_NL + norm(dv);

    % Apply impulse
    v_curr_NL = v_curr_NL + dv;

    % --- Nonlinear propagation ---
    s_temp = [r_curr_NL; v_curr_NL];

    [~, s_out] = ode45(f_nonlin, [0 dt_imp], s_temp, opts);

    s_next = s_out(end,:)';

    r_curr_NL = s_next(1:3);
    v_curr_NL = s_next(4:6);

    traj_vbar_NL(k,:) = r_curr_NL';
end

fprintf('--- POINT 5.4 NONLINEAR CHECK ---\n');
fprintf('Total Delta-V (Nonlinear) = %.6f m/s\n\n', DeltaV_vbar_NL*1000);

% =========================================================
%  ERROR HCW vs NONLINEAR (trajectory)
% ==========================================================

err_vbar = vecnorm(traj_vbar_NL - traj_vbar,2,2);

fprintf('Final position error (HCW) = %.6f m\n', norm(traj_vbar(end,:))*1000);
fprintf('Final position error (NL)  = %.6f m\n', norm(traj_vbar_NL(end,:))*1000);

% ======= TRAJECTORY PLOT =====================================

figure('Name','V-bar HCW vs Nonlinear','Color','w');
hold on; grid on; axis equal;

plot3(traj_vbar_NL(:,1)*1000, traj_vbar_NL(:,2)*1000, traj_vbar_NL(:,3)*1000, ...
    'r-o','LineWidth',1.5,'DisplayName','Nonlinear');

plot3(traj_vbar(:,1)*1000, traj_vbar(:,2)*1000, traj_vbar(:,3)*1000, ...
    'b--o','LineWidth',1,'DisplayName','Trajectory');

plot3(0,0,0,'ks','MarkerFaceColor','k','DisplayName','Docking Port');

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('V-bar Approach: HCW vs Nonlinear');
legend('Location','best');
view(40,25);

% ======= ERROR PLOT =====================================

figure('Name','Error V-bar HCW vs Nonlinear','Color','w');
plot(1:N_imp, err_vbar*1000,'k','LineWidth',1.5);
grid on;

xlabel('Impulse step');
ylabel('Error [m]');
title('Deviation HCW vs Nonlinear – V-bar Approach');


% Function to evaluate the final position mismatch w. r. t. the target
% after nonlinear propagation for a given first impulsive maneuver d_V1

function err = final_position_error(dv1, r0, v0, t_man, f_nonlin, opts, r_target)
s0 = [r0; v0 + dv1];
[~, s] = ode45(f_nonlin, [0 t_man], s0, opts);
r_final = s(end,1:3)';
err = r_final - r_target;
end