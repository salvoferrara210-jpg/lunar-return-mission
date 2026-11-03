clear all;
close all;
clc;

%% =========================================================
%  PART 3 – Perturbed Orbital Dynamics
%
% To assess the impact of the lunar orbital environment
% on the LM trajectory and rendezvous.
%
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
J2 = 2.033e-4;

v_LM = sqrt(mu_M / r_LM);   % Circular velocity on LM orbit [km/s]
v_MS = sqrt(mu_M / r_MS);   % Circular velocity on MS orbit [km/s]

a = (r_LM + r_MS) / 2;                   % Transfer ellipse SMA [km]
v_perigee = sqrt(mu_M * (2/r_LM - 1/a)); % Velocity at perigee [km/s]
v_apogee  = sqrt(mu_M * (2/r_MS - 1/a)); % Velocity at apogee  [km/s]

delta_V1      = v_perigee - v_LM;    % First impulse  [km/s]
delta_V2      = v_MS - v_apogee;     % Second impulse [km/s]
total_delta_V = delta_V1 + delta_V2; % Total delta-V  [km/s]
delta_t = pi * sqrt(a^3 / mu_M);     % Transfer time of flight [s]

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
T_tot = t_wait + delta_t;

%% =========================================================
%  POINT 2 – Implement the J2 perturbation 
% ==========================================================

% LM initial state
r_LM0 = [r_LM; 0; 0];
v_LM0 = [0; v_LM; 0];

% MS initial state 
r_MS0 = [r_MS*cos(phi0); r_MS*sin(phi0); 0];
v_MS0 = [-v_MS*sin(phi0); v_MS*cos(phi0); 0];

% State vector initial conditions t = 0
x_LM0 = [r_LM0;v_LM0];
x_MS0= [r_MS0; v_MS0];

% --- Two Body Problem ODE ---
two_body = @(~, s) [s(4);
                    s(5);
                    s(6);
                   -mu_M * s(1) / (s(1)^2 + s(2)^2 + s(3)^2)^1.5;
                   -mu_M * s(2) / (s(1)^2 + s(2)^2 + s(3)^2)^1.5;
                   -mu_M * s(3) / (s(1)^2 + s(2)^2 + s(3)^2)^1.5];

% --- Cowell Method ODE ---
cowell = @(~, s) [s(4:6);
                 -mu_M / norm(s(1:3))^3 * s(1:3) + ...
                  3/2 * J2 * mu_M * R_M^2 / norm(s(1:3))^5 * [s(1)*(5*(s(3)/norm(s(1:3)))^2 - 1); ...
                                                              s(2)*(5*(s(3)/norm(s(1:3)))^2 - 1); ...
                                                              s(3)*(5*(s(3)/norm(s(1:3)))^2 - 3)] ];

% Options
opts   = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% Phase 1 : Waiting Time on Initial Orbits
N_wait = 6000;
tspan_wait = linspace(0,t_wait,N_wait);

% Unperturbed
[t_kep_w, X_LM_kep_w] = ode45(two_body,tspan_wait,x_LM0,opts);
[~,     X_MS_kep_w] = ode45(two_body,tspan_wait,x_MS0,opts);

% Perturbed
[t_J2_w, X_LM_J2_w] = ode45(cowell,tspan_wait,x_LM0,opts);
[~,     X_MS_J2_w] = ode45(cowell,tspan_wait,x_MS0,opts);

% ==========================================================
%  IDEAL KEPLERIAN MANEUVER
% ==========================================================

% Calculate the required velocity at the perigee 
v_trans = sqrt(mu_M*(2/r_LM - 1/a));

% Position vector of LM at the departure point
r_LM_dep = X_LM_kep_w(end,1:3)';

% Local tangential unit vector
t_hat = [-r_LM_dep(2); r_LM_dep(1); 0] / norm(r_LM_dep(1:2)); % tangential direction

% Delta_V vector 
dV1_vec = (v_trans - v_LM)*t_hat;

% Intial State for transfer phase
X_LM_t0_kep = [X_LM_kep_w(end,1:3)'; X_LM_kep_w(end,4:6)'+dV1_vec];
X_MS_t0_kep = X_MS_kep_w(end,:)';

% ==========================================================
% J2 PERTURBATED MANEUVER
% ==========================================================

% Actual velocity vector of the perturbed LM at the departure point
v_J2_curr = X_LM_J2_w(end,4:6)';

% True fligth path direction
t_hat_J2 = v_J2_curr / norm(v_J2_curr);

% Delta_V vector
dV1_vec_J2 = (v_trans - v_LM) * t_hat_J2;

% Initial State for the pertubed transfer phase
X_LM_t0_J2 = [X_LM_J2_w(end,1:3)'; X_LM_J2_w(end,4:6)' + dV1_vec_J2];
X_MS_t0_J2 = X_MS_J2_w(end,:)';

% Phase 2: Hohmann Transfer
N_trans = 6000;
tspan_trans = linspace(t_wait, t_wait+delta_t,N_trans);

% Unperturbed
[t_tr_kep, X_LM_tr_kep] = ode45(two_body,tspan_trans,X_LM_t0_kep,opts);
[~, X_MS_tr_kep] = ode45(two_body,tspan_trans,X_MS_t0_kep,opts);

% Perturbed
[t_tr_J2, X_LM_tr_J2] = ode45(cowell,tspan_trans,X_LM_t0_J2,opts);
[~, X_MS_tr_J2] = ode45(cowell,tspan_trans,X_MS_t0_J2,opts);

% --- Concatenate the Trajectories ---
% Unperturbed
t_kep_full = [t_kep_w; t_tr_kep(2:end)];
X_LM_kep_full = [X_LM_kep_w; X_LM_tr_kep(2:end,:)];
X_MS_kep_full = [X_MS_kep_w; X_MS_tr_kep(2:end,:)];

% Perturbed
t_J2_full = [t_J2_w; t_tr_J2(2:end)];
X_LM_J2_full = [X_LM_J2_w; X_LM_tr_J2(2:end,:)];
X_MS_J2_full = [X_MS_J2_w; X_MS_tr_J2(2:end,:)];

nw_kep = length(t_kep_w); % Number of the point in the ideal Keplerian Model
nw_J2 = length(t_J2_w);   % Number of the point in the Cowell Model

% Compute orbital element along trajectories
OE_LM_kep = states2OE(X_LM_kep_full, mu_M);
OE_LM_j2  = states2OE(X_LM_J2_full,  mu_M);
OE_MS_kep = states2OE(X_MS_kep_full, mu_M);
OE_MS_j2  = states2OE(X_MS_J2_full,  mu_M);

% Convert time
t_h_kep = t_kep_full / 3600;
t_h_J2 = t_J2_full / 3600;
t_wait_h = t_wait / 3600;

%% =========================================================================
%  FIGURE 1 – LM Trajectory
%  =========================================================================

% Plot Moon
theta_c = linspace(0, 2*pi, 300);
R_M7 = R_M*0.7;
x_moon = R_M7*cos(theta_c); 
y_moon = R_M7*sin(theta_c);
x_c1   = r_LM*cos(theta_c);  
y_c1   = r_LM*sin(theta_c);
x_c2   = r_MS*cos(theta_c);  
y_c2   = r_MS*sin(theta_c);
 
figure('Name','Point 3.2 – LM Trajectory','Position',[50 50 900 750]);
hold on; axis equal; grid on;
fill(x_moon, y_moon, [0.75 0.75 0.75], 'EdgeColor','none', 'DisplayName','Moon');
plot(x_c1, y_c1, 'b:', 'LineWidth',0.8, 'HandleVisibility','off');
plot(x_c2, y_c2, 'k:', 'LineWidth',0.8, 'HandleVisibility','off');
 
% Waiting arc
plot(X_LM_kep_full(1:nw_kep,1), X_LM_kep_full(1:nw_kep,2), ...
     'Color',[0 0.6 0], 'LineWidth',1.6, 'DisplayName','LM wait – Keplerian');
plot(X_LM_J2_full(1:nw_J2,1),   X_LM_J2_full(1:nw_J2,2), ...
     '--', 'Color',[0.4 0.8 0.2], 'LineWidth',1.6, 'DisplayName','LM wait – J2');
 
% Transfer arc
plot(X_LM_kep_full(nw_kep:end,1), X_LM_kep_full(nw_kep:end,2), ...
     'b-', 'LineWidth',1.8, 'DisplayName','LM transfer – Keplerian');
plot(X_LM_J2_full(nw_J2:end,1),  X_LM_J2_full(nw_J2:end,2), ...
     'r--', 'LineWidth',1.8, 'DisplayName','LM transfer – J2');
 
% Markers
plot(X_LM_kep_full(1,1),   X_LM_kep_full(1,2),   'gs', 'MarkerSize',10, 'MarkerFaceColor','g', 'DisplayName','LM t=0');
plot(X_LM_kep_full(nw_kep,1), X_LM_kep_full(nw_kep,2), 'k^', 'MarkerSize',10, 'MarkerFaceColor','y', 'DisplayName','LM at \DeltaV_1');
plot(X_LM_kep_full(end,1), X_LM_kep_full(end,2), 'b^', 'MarkerSize',10, 'MarkerFaceColor','c', 'DisplayName','LM arrival (Kep)');
plot(X_LM_J2_full(end,1),  X_LM_J2_full(end,2),  'r^', 'MarkerSize',10, 'MarkerFaceColor','r', 'DisplayName','LM arrival (J2)');
 
xlabel('X [km]','FontSize',12); ylabel('Y [km]','FontSize',12);
title('Point 3.2 – LM Trajectory: Keplerian vs J2 (wait + transfer)','FontSize',13);
legend('Location','best','FontSize',9); set(gca,'FontSize',11);
 
%% =========================================================================
%  FIGURE 2 – MS Trajectory
%  =========================================================================
figure('Name','Point 3.2 – MS Trajectory','Position',[100 80 900 750]);
hold on; axis equal; grid on;
fill(x_moon, y_moon, [0.75 0.75 0.75], 'EdgeColor','none', 'DisplayName','Moon');
plot(x_c2, y_c2, 'k:', 'LineWidth',0.8, 'HandleVisibility','off');
 
plot(X_MS_kep_full(:,1), X_MS_kep_full(:,2), 'b-',  'LineWidth',1.8, 'DisplayName','MS – Keplerian');
plot(X_MS_J2_full(:,1),  X_MS_J2_full(:,2),  'r--', 'LineWidth',1.8, 'DisplayName','MS – J2');
 
plot(X_MS_kep_full(1,1),   X_MS_kep_full(1,2),   'bs', 'MarkerSize',10, 'MarkerFaceColor','b', 'DisplayName','MS t=0');
plot(X_MS_kep_full(end,1), X_MS_kep_full(end,2), 'b^', 'MarkerSize',10, 'MarkerFaceColor','c', 'DisplayName','MS end (Kep)');
plot(X_MS_J2_full(end,1),  X_MS_J2_full(end,2),  'r^', 'MarkerSize',10, 'MarkerFaceColor','r', 'DisplayName','MS end (J2)');
 
xlabel('X [km]','FontSize',12); ylabel('Y [km]','FontSize',12);
title('Point 3.2 – MS Trajectory: Keplerian vs J2 (full duration)','FontSize',13);
legend('Location','best','FontSize',9); set(gca,'FontSize',11);
 
%% =========================================================================
%  FIGURE 3 – Orbital Elements: LM
%  =========================================================================
OE_labels = {'a  [km]','e  [-]','i  [deg]','\Omega  [deg]','\omega  [deg]','\nu  [deg]'};
OE_titles = {'Semi-major axis','Eccentricity','Inclination','RAAN','Arg. of perigee','True anomaly'};

% States accurate at key moments
idx_wait = nw_kep;
idx_end = size(X_LM_kep_full,1);

% At the end of the wait
OE_LM_kep_wait = states2OE(X_LM_kep_full(idx_wait,:), mu_M);
OE_LM_J2_wait = states2OE(X_LM_J2_full(idx_wait,:), mu_M);
dOE_LM_wait = OE_LM_J2_wait - OE_LM_kep_wait;

% At the final Rendezvous
OE_LM_kep_end = states2OE(X_LM_kep_full(idx_end,:), mu_M);
OE_LM_J2_end = states2OE(X_LM_J2_full(idx_end,:), mu_M);
dOE_LM_end = OE_LM_J2_end - OE_LM_kep_end;

figure('Name','Point 3.2 – OE LM','Position',[150 50 1200 800]);
for k = 1:6
    subplot(3,2,k); hold on; grid on;
    xline(t_wait_h,'k--','LineWidth',1.2,'HandleVisibility','off');
    plot(t_h_kep, OE_LM_kep(:,k), 'b-',  'LineWidth',1.5, 'DisplayName','Keplerian');
    plot(t_h_J2,  OE_LM_j2(:,k),  'r--', 'LineWidth',1.5, 'DisplayName','J2 perturbed');
    xlabel('Time [h]','FontSize',10);
    ylabel(OE_labels{k},'FontSize',10);
    title(['LM – ', OE_titles{k}],'FontSize',10);
    legend('Location','best','FontSize',8); set(gca,'FontSize',9);
end
sgtitle('Point 3.2 – Orbital Elements: LM   (dashed line = \DeltaV_1)', ...
        'FontSize',13,'FontWeight','bold');
 
%% =========================================================================
%  FIGURE 4 – Orbital Elements: MS
%  =========================================================================

% At the end of the wait
OE_MS_kep_wait = states2OE(X_MS_kep_full(idx_wait,:), mu_M);
OE_MS_J2_wait = states2OE(X_MS_J2_full(idx_wait,:), mu_M);
dOE_MS_wait = OE_MS_J2_wait - OE_MS_kep_wait;

% At the final Rendezvous
OE_MS_kep_end = states2OE(X_MS_kep_full(idx_end,:), mu_M);
OE_MS_J2_end = states2OE(X_MS_J2_full(idx_end,:), mu_M);
dOE_MS_end = OE_MS_J2_end - OE_MS_kep_end;

figure('Name','Point 3.2 – OE MS','Position',[200 80 1200 800]);
for k = 1:6
    subplot(3,2,k); hold on; grid on;
    xline(t_wait_h,'k--','LineWidth',1.2,'HandleVisibility','off');
    plot(t_h_kep, OE_MS_kep(:,k), 'b-',  'LineWidth',1.5, 'DisplayName','Keplerian');
    plot(t_h_J2,  OE_MS_j2(:,k),  'r--', 'LineWidth',1.5, 'DisplayName','J2 perturbed');
    xlabel('Time [h]','FontSize',10);
    ylabel(OE_labels{k},'FontSize',10);
    title(['MS – ', OE_titles{k}],'FontSize',10);
    legend('Location','best','FontSize',8); set(gca,'FontSize',9);
end
sgtitle('Point 3.2 – Orbital Elements: MS   (dashed line = \DeltaV_1)', ...
        'FontSize',13,'FontWeight','bold');

fprintf('\n========================================================================================\n');
fprintf(' PART 3.2: ORBITAL ELEMENTS DRIFT SUMMARY (Perturbed J2 - Ideal Keplerian)\n');
fprintf('========================================================================================\n');

fprintf('LUNAR MODULE (LM):\n');
fprintf('                 Delta a [km]   Delta e [-]     Delta i [deg]   Delta RAAN [deg]  Delta w [deg]   Delta nu [deg]\n');
fprintf(' - End of Wait: %12.6f %14.6e %14.6f %16.6f %16.6f %16.6f\n', dOE_LM_wait);
fprintf(' - Rendezvous : %12.6f %14.6e %14.6f %16.6f %16.6f %16.6f\n\n', dOE_LM_end);

fprintf('MOTHERSHIP (MS):\n');
fprintf('                 Delta a [km]   Delta e [-]     Delta i [deg]   Delta RAAN [deg]  Delta w [deg]   Delta nu [deg]\n');
fprintf(' - End of Wait: %12.6f %14.6e %14.6f %16.6f %16.6f %16.6f\n', dOE_MS_wait);
fprintf(' - Rendezvous : %12.6f %14.6e %14.6f %16.6f %16.6f %16.6f\n', dOE_MS_end);
fprintf('========================================================================================\n\n');

%% =========================================================================
%  FIGURE 5 – Position drift due to J2
%  =========================================================================
t_common   = t_J2_full;
X_LM_kep_i = interp1(t_kep_full, X_LM_kep_full, t_common, 'pchip');
X_MS_kep_i = interp1(t_kep_full, X_MS_kep_full, t_common, 'pchip');
 
dr_LM = vecnorm(X_LM_J2_full(:,1:3) - X_LM_kep_i(:,1:3), 2, 2) * 1e3;  % [m]
dr_MS = vecnorm(X_MS_J2_full(:,1:3) - X_MS_kep_i(:,1:3), 2, 2) * 1e3;  % [m]

idx_wait = nw_kep;
idx_end = length(dr_LM);

err_LM_wait = dr_LM(idx_wait);
err_LM_end = dr_LM(idx_end);
max_err_LM = max(dr_LM);

err_MS_wait = dr_MS(idx_wait);
err_MS_end = dr_MS(idx_end);
max_err_MS = max(dr_MS);

fprintf('\n======================================================\n');
fprintf(' PART 3.2: J2 POSITION DRIFT SUMMARY (Perturbed vs Kep)\n');
fprintf('======================================================\n');

fprintf('LUNAR MODULE (LM):\n');
fprintf(' - Drift at end of Wait Phase: %8.2f m\n', err_LM_wait);
fprintf(' - Drift at Final Rendezvous:  %8.2f m\n', err_LM_end);
fprintf(' - Maximum Drift overall:      %8.2f m\n\n', max_err_LM);

fprintf('MOTHERSHIP (MS):\n');
fprintf(' - Drift at end of Wait Phase: %8.2f m\n', err_MS_wait);
fprintf(' - Drift at Final Rendezvous:  %8.2f m\n', err_MS_end);
fprintf(' - Maximum Drift overall:      %8.2f m\n', max_err_MS);

fprintf('======================================================\n\n');
figure('Name','Point 3.2 – J2 Position Drift','Position',[250 100 1000 430]);
subplot(1,2,1); hold on; grid on;
xline(t_wait_h,'k--','LineWidth',1.2,'HandleVisibility','off');
plot(t_common/3600, dr_LM, 'b-', 'LineWidth',1.8);
xlabel('Time [h]','FontSize',12); ylabel('||\Deltar|| [m]','FontSize',12);
title('J2 Position Drift – LM','FontSize',12); set(gca,'FontSize',11);
 
subplot(1,2,2); hold on; grid on;
xline(t_wait_h,'k--','LineWidth',1.2,'HandleVisibility','off');
plot(t_common/3600, dr_MS, 'r-', 'LineWidth',1.8);
xlabel('Time [h]','FontSize',12); ylabel('||\Deltar|| [m]','FontSize',12);
title('J2 Position Drift – MS','FontSize',12); set(gca,'FontSize',11);
sgtitle('Point 3.2 – Position Drift due to J2  (dashed = \DeltaV_1)', ...
        'FontSize',13,'FontWeight','bold');
 
%% =========================================================================
%  FUNCTION - Calculation of Orbital Elements (OE)
%  =========================================================================

function OE = states2OE(X, mu_M)
    N = size(X,1); 
    OE = zeros(N,6);

    for k = 1:N
        r_vec = X(k,1:3)'; % Position vector
        v_vec = X(k,4:6)'; % Vecocity vector
        r = norm(r_vec); 
        v = norm(v_vec);
        a     = -mu_M / (2*(v^2/2 - mu_M/r)); % Semi-major axis
        h_vec = cross(r_vec, v_vec);          % Angular momentum vector
        h = norm(h_vec);
        e_vec = (1/mu_M)*((v^2 - mu_M/r)*r_vec - dot(r_vec,v_vec)*v_vec); % Eccentricity vector
        e     = norm(e_vec);
        i     = acosd(h_vec(3)/h);     % Orbital plane inclination
        N_vec = cross([0;0;1], h_vec); % Axis of nodes vector
        N_n = norm(N_vec);

        if N_n < 1e-10, RAAN = 0;
        else, RAAN = acosd(N_vec(1)/N_n); % Right ascension of the ascending node
            if N_vec(2)<0, RAAN=360-RAAN; 
            end; end
        if N_n<1e-10 || e<1e-10, omega = 0;
        else, omega=acosd(max(-1,min(1,dot(N_vec,e_vec)/(N_n*e)))); % Argument of periapsis
            if e_vec(3)<0, omega=360-omega;
            end; end
        if e<1e-10, nu=0;
        else, nu=acosd(max(-1,min(1,dot(e_vec,r_vec)/(e*r)))); % True anomaly 
            if dot(r_vec,v_vec)<0, nu=360-nu;
            end; end
        OE(k,:) = [a, e, i, RAAN, omega, nu];
    end
end

%% =========================================================
%  POINT 3 – Third-body perturbation
%  =========================================================

mu_E = 3.986e5; % [km^3/s^2]
d_EM = 384400; % [km]
n_E = sqrt((mu_E + mu_M) / d_EM^3); % Mean Motion

cowell_EM = @(t,s) cowell_E_M(t, s, mu_M, R_M, J2, mu_E, d_EM, n_E);

% Phase 1: Waiting Time
[~, X_LM_wait_EM] = ode45(cowell_EM,tspan_wait,x_LM0,opts);
[~, X_MS_wait_EM] = ode45(cowell_EM,tspan_wait,x_MS0,opts);

% Phase 2: Trasfer
v_curr_EM = X_LM_wait_EM(end,4:6)';
t_hat_EM = v_curr_EM / norm(v_curr_EM);
dV1_vec_EM = (v_trans - v_LM)* t_hat_EM;

X_LM_t0_EM = [X_LM_wait_EM(end,1:3)'; X_LM_wait_EM(end,4:6)' + dV1_vec_EM];
X_MS_t0_EM = [X_MS_wait_EM(end,:)'];

[~, X_LM_tr_EM] = ode45(cowell_EM,tspan_trans,X_LM_t0_EM,opts);
[~, X_MS_tr_EM] = ode45(cowell_EM,tspan_trans,X_MS_t0_EM,opts);

% --- Concatenate the Trajectories ---
X_LM_EM_full = [X_LM_wait_EM; X_LM_tr_EM(2:end,:)];
X_MS_EM_full = [X_MS_wait_EM; X_MS_tr_EM(2:end,:)];

% Ideal relative position vector at the time of Rendezvous
rel_pos_kep = X_LM_kep_full(end,1:3) - X_MS_kep_full(end,1:3);

% Real relative position vector (Perturbated J2+3B)
rel_pos_EM = X_LM_EM_full(end,1:3) - X_MS_EM_full(end,1:3);

% Relative drift
relative_drift = norm(rel_pos_EM - rel_pos_kep); % [km]
rel_drift_m = relative_drift * 1e3; % [m]

fprintf('\n======================================================\n');
fprintf(' PART 3.3 \n');
fprintf('======================================================\n');
fprintf('Relative Position Drift (LM vs MS) with J2 + 3B: %.2f km\n\n', relative_drift);

%% =========================================================================
%  FIGURE 6 – Point 3.3: Trajectories J2+3B vs Keplerian
%  =========================================================================
% Interpolate Keplerian solution onto the J2+3B time grid
% (J2+3B and Keplerian grids may differ due to adaptive stepping)
t_EM_full = [tspan_wait(:); tspan_trans(2:end)'];

X_LM_kep_i = interp1(t_kep_full, X_LM_kep_full, t_EM_full, 'pchip');
X_MS_kep_i = interp1(t_kep_full, X_MS_kep_full, t_EM_full, 'pchip');

% Absolute position drift of LM and MS under J2+3B vs Keplerian
dr_LM_EM = vecnorm(X_LM_EM_full(:,1:3) - X_LM_kep_i(:,1:3), 2, 2) * 1e3; % [m]
dr_MS_EM = vecnorm(X_MS_EM_full(:,1:3) - X_MS_kep_i(:,1:3), 2, 2) * 1e3; % [m]

figure('Name','Point 3.3 – J2+3B Position Drift','Position',[350 80 700 430]);
hold on; grid on;
xline(t_wait_h, 'k--', 'LineWidth',1.2, 'DisplayName','\DeltaV_1 burn');
plot(t_EM_full/3600, dr_LM_EM, 'b-',  'LineWidth',1.8, 'DisplayName','LM drift');
plot(t_EM_full/3600, dr_MS_EM, 'r--', 'LineWidth',1.8, 'DisplayName','MS drift');
xlabel('Time [h]','FontSize',12);
ylabel('||\Deltar|| [m]','FontSize',12);
title('Point 3.3 – Position Drift: J2+3B vs Keplerian','FontSize',13);
legend('Location','best','FontSize',9);

%  ==================================================================
%  FUNCTION - Cowell Method with J2 and Earth Third-Body Perturbation
% ===================================================================

function ds = cowell_E_M(t, s, mu_M, R_M, J2, mu_E, d_EM, n_E)
r_vec = s(1:3);
v_vec = s(4:6);
r = norm(r_vec);

a_kep = -mu_M / r^3 * r_vec; % Keplerian acceleration

fac = (3/2) * J2 * (mu_M * R_M^2) / r^5;
a_J2_x = fac * r_vec(1) * (5 * (r_vec(3)/r)^2 - 1);
a_J2_y = fac * r_vec(2) * (5 * (r_vec(3)/r)^2 - 1);
a_J2_z = fac * r_vec(3) * (5 * (r_vec(3)/r)^2 - 3);

a_J2 = [a_J2_x;a_J2_y;a_J2_z]; % J2 acceleration

% Position Vector Earth
r_E = d_EM * [cos(n_E*t); sin(n_E*t);0];
% Position from satellite to Earth vector
r_st_E = r_E - r_vec;

% Perturbed acceleration of the Third Body
a_3B = mu_E * ((r_st_E / norm(r_st_E)^3) - (r_E / norm(r_E)^3));

ds = [v_vec; a_kep + a_J2 + a_3B];
end