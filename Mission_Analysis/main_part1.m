clear all;
close all;
clc;

%% =========================================================
%  PART 1 – Orbital Transfer 
%
% Analytically compute the minium delta-V to transfer the 
% Lunar Module (LM) from its initial Low Lunar Orbit (LLO) 
% to the mothership's (MS) orbit, accounting for the initial 
% phasing angle (phi_0)
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

%% =========================================================
%  POINT 1 – Hohmann Transfer: Delta-V and Time of Flight
% ==========================================================
v_LM = sqrt(mu_M / r_LM);   % Circular velocity on LM orbit [km/s]
v_MS = sqrt(mu_M / r_MS);   % Circular velocity on MS orbit [km/s]

a         = (r_LM + r_MS) / 2;                   % Transfer ellipse SMA [km]
v_perigee = sqrt(mu_M * (2/r_LM - 1/a));         % Velocity at perigee [km/s]
v_apogee  = sqrt(mu_M * (2/r_MS - 1/a));         % Velocity at apogee  [km/s]

delta_V1      = v_perigee - v_LM;                % First impulse  [km/s]
delta_V2      = v_MS - v_apogee;                 % Second impulse [km/s]
total_delta_V = delta_V1 + delta_V2;             % Total delta-V  [km/s]
delta_t = pi * sqrt(a^3 / mu_M);                 % Transfer time of flight [s]

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
title('Point 1.3 – Hohmann Transfer Geometry','FontSize',13);
legend('Location','northeastoutside','FontSize',8);
xlim([-(r_MS+500),(r_MS+500)]); ylim([-(r_MS+500),(r_MS+500)]);
hold off;

%% =========================================================
%  POINT 4 – Numerical Propagation & Verification (Full Maneuver)
% ==========================================================
% --- Two-body ODE ---
two_body = @(~, s) [s(3);
                    s(4);
                   -mu_M * s(1) / (s(1)^2 + s(2)^2)^1.5;
                   -mu_M * s(2) / (s(1)^2 + s(2)^2)^1.5];

opts   = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% --- PHASE 1: Waiting Time on Initial Orbits ---
% Initial conditions at t=0
s0_LM_wait = [r_LM; 0; 0; v_LM];
s0_MS_wait = [r_MS * cos(phi0); r_MS * sin(phi0); -v_MS * sin(phi0); v_MS * cos(phi0)];

N_wait = 6000;
t_span_wait = linspace(0, t_wait, N_wait);
[t_wait_num, sol_LM_wait] = ode45(two_body, t_span_wait, s0_LM_wait, opts);
[~,          sol_MS_wait] = ode45(two_body, t_span_wait, s0_MS_wait, opts);

v_dir     = s_LM_burn(3:4) / norm(s_LM_burn(3:4));
s0_LM_trans = [s_LM_burn(1:2); s_LM_burn(3:4) + delta_V1 * v_dir];

% MS continues on its circular orbit (no delta-V)
s0_MS_trans = sol_MS_wait(end, :)';

N_trans = 6000;
t_span_trans = linspace(t_wait, t_wait + delta_t, N_trans);
[t_trans_num, sol_LM_trans] = ode45(two_body, t_span_trans, s0_LM_trans, opts);
[~,           sol_MS_trans] = ode45(two_body, t_span_trans, s0_MS_trans, opts);

% --- Concatenate the Trajectories ---
t_num  = [t_wait_num; t_trans_num(2:end)];
sol_LM = [sol_LM_wait; sol_LM_trans(2:end, :)];
sol_MS = [sol_MS_wait; sol_MS_trans(2:end, :)];
N_tot  = length(t_num);

% --- Analytical Solution for the Full Maneuver ---
x_an_LM = zeros(N_tot, 1); y_an_LM = zeros(N_tot, 1);
x_an_MS = zeros(N_tot, 1); y_an_MS = zeros(N_tot, 1);

n_tr       = sqrt(mu_M / a^3);
theta_burn = n_LM * t_wait; 

for k = 1:N_tot
    t = t_num(k);
    
    % MS is always in circular orbit
    theta_ms = phi0 + n_MS * t;
    x_an_MS(k) = r_MS * cos(theta_ms);
    y_an_MS(k) = r_MS * sin(theta_ms);
    
    if t <= t_wait
        % Phase 1: LM is in circular orbit
        theta_lm = n_LM * t;
        x_an_LM(k) = r_LM * cos(theta_lm);
        y_an_LM(k) = r_LM * sin(theta_lm);
    else
        % Phase 2: LM is on the transfer ellipse
        dt_trans = t - t_wait;
        M = n_tr * dt_trans;
        Ek = M;
        % Newton-Raphson for Kepler's Equation
        for iter = 1:100
            f_eq  = Ek - ecc*sin(Ek) - M;
            df_eq = 1  - ecc*cos(Ek);
            Ek = Ek - f_eq/df_eq;
            if abs(f_eq/df_eq) < 1e-14, break; end
        end
        % True anomaly f (since E is 0 to pi, f is 0 to pi)
        f_true = 2 * atan( sqrt((1+ecc)/(1-ecc)) * tan(Ek/2) );
        r_mag  = a * (1 - ecc * cos(Ek));
        
        % Total inertial angle
        theta_lm = theta_burn + f_true;
        x_an_LM(k) = r_mag * cos(theta_lm);
        y_an_LM(k) = r_mag * sin(theta_lm);
    end
end

% --- Absolute position errors ---
err_LM = sqrt((sol_LM(:,1)-x_an_LM).^2 + (sol_LM(:,2)-y_an_LM).^2);
err_MS = sqrt((sol_MS(:,1)-x_an_MS).^2 + (sol_MS(:,2)-y_an_MS).^2);

% ---- Relative position errors ---
delta_r_num_x = sol_LM(:,1) - sol_MS(:,1);
delta_r_num_y = sol_LM(:,2) - sol_MS(:,2);
delta_r_an_x  = x_an_LM    - x_an_MS;
delta_r_an_y  = y_an_LM    - y_an_MS;
err_rel = sqrt((delta_r_num_x - delta_r_an_x).^2 + (delta_r_num_y - delta_r_an_y).^2);   % [km]

% --- Rendezvous final position check ---
r_LM_final = sol_LM(end, 1:2);
r_MS_final = sol_MS(end, 1:2);
pos_error_final = norm(r_LM_final - r_MS_final);   % [km]

fprintf('----- RENDEZVOUS CHECK -----\n');
fprintf('Final position error  = %.6e km\n', pos_error_final);
tol_pos = 1e-6;
if pos_error_final < tol_pos
    fprintf('Position match       --> PASS\n');
else
    fprintf('Position match       --> FAIL\n');
end

% --- Explicit verification against tolerance ---
tol = 1e-6;   % [km] = 1 mm  
fprintf('===== POINT 4 =====\n');
fprintf('Integration tolerances : RelTol = 1e-12,  AbsTol = 1e-14\n');
fprintf('Verification tolerance : %.0e km  (sub-millimetre)\n', tol);
fprintf('Max |err_LM|           : %.4e km  -->  %s\n', max(err_LM),  verdict(max(err_LM),  tol));
fprintf('Max |err_MS|           : %.4e km  -->  %s\n', max(err_MS),  verdict(max(err_MS),  tol));
fprintf('Max |err_rel|          : %.4e km  -->  %s\n\n', max(err_rel), verdict(max(err_rel), tol));

t_h = t_num / 3600;   % time axis in hours
t_wait_h = t_wait / 3600;

% --- Figure 1: error plots ---
figure('Name','Point 1.4 – Position Errors','Color','w');
subplot(2,1,1);
semilogy(t_h, err_LM, 'b-',  'LineWidth',1.3,'DisplayName','LM  ||r_{num} - r_{an}||');
hold on;
semilogy(t_h, err_MS, 'r--', 'LineWidth',1.3,'DisplayName','MS  ||r_{num} - r_{an}||');
xline(t_wait_h, 'm-.', 'LineWidth',1.0, 'HandleVisibility','off');
yline(tol,'k:','LineWidth',1.2, ...
    'Label',sprintf('tol = %.0e km', tol), ...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
grid on;
xlabel('Time  [h]','FontWeight','bold');
ylabel('Absolute error  [km]','FontWeight','bold');
title('Absolute Position Error – LM and MS');
legend('Location','best','FontSize',9);
hold off;

subplot(2,1,2);
semilogy(t_h, err_rel, 'k-', 'LineWidth',1.3, ...
    'DisplayName','||\Deltar_{num}(t) - \Deltar_{an}(t)||');
hold on;
xline(t_wait_h, 'm-.', 'LineWidth',1.0, 'DisplayName', '\DeltaV_1 Burn Time');
yline(tol,'k:','LineWidth',1.2, ...
    'Label',sprintf('tol = %.0e km', tol), ...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
grid on;
xlabel('Time  [h]','FontWeight','bold');
ylabel('Relative position error  [km]','FontWeight','bold');
title('Relative Position Error  ||\Deltar_{num}(t) - \Deltar_{an}(t)||');
legend('Location','best','FontSize',9);
hold off;
sgtitle('Point 1.4 – Propagation Verification (Full Maneuver)','FontSize',13);


%% =========================================================
%  POINT 5 – Kepler Constants of Motion for LM
% ==========================================================
eps_num = zeros(N_tot,1);
h_num   = zeros(N_tot,1);
ecc_num = zeros(N_tot,1);

eps_ref = zeros(N_tot,1);
h_ref   = zeros(N_tot,1);
ecc_ref = zeros(N_tot,1);

for k = 1:N_tot
    rx = sol_LM(k,1);  
    ry = sol_LM(k,2);
    vx = sol_LM(k,3);  
    vy = sol_LM(k,4);
    
    r_mag = sqrt(rx^2 + ry^2);
    v_mag = sqrt(vx^2 + vy^2);
    
    % Numerical constants
    eps_num(k) = 0.5*v_mag^2 - mu_M/r_mag;
    hz         = rx*vy - ry*vx;
    h_num(k)   = abs(hz);
    ex         = ( vy*hz)/mu_M - rx/r_mag;
    ey         = (-vx*hz)/mu_M - ry/r_mag;
    ecc_num(k) = sqrt(ex^2 + ey^2);
    
    % Analytical references 
    if t_num(k) <= t_wait
        eps_ref(k) = -mu_M / (2*r_LM);
        h_ref(k)   = r_LM * v_LM;
        ecc_ref(k) = 0;
    else
        eps_ref(k) = -mu_M / (2*a);
        h_ref(k)   = sqrt(mu_M * p);
        ecc_ref(k) = ecc;
    end
end

fprintf('===== POINT 5 =====\n');
fprintf('Max numerical drift from analytical constant (over full maneuver):\n');
fprintf('  delta_eps : %.4e km^2/s^2\n', max(abs(eps_num - eps_ref)));
fprintf('  delta_h   : %.4e km^2/s\n',   max(abs(h_num   - h_ref)));
fprintf('  delta_ecc : %.4e\n',           max(abs(ecc_num - ecc_ref)));

figure('Name','Point 1.5 – Kepler Constants of Motion','Color','w');
subplot(3,1,1);
plot(t_h, eps_num, 'b-', 'LineWidth',2.0,'DisplayName','Numerical \epsilon');
hold on;
plot(t_h, eps_ref, 'r--', 'LineWidth',1.2,'DisplayName','Analytic Ref.');
xline(t_wait_h, 'k:', 'LineWidth',1.0, 'HandleVisibility','off');
grid on;
xlabel('Time [h]','FontWeight','bold'); ylabel('\epsilon  [km^2/s^2]','FontWeight','bold');
title('Specific Orbital Energy  \epsilon');
legend('Location','best');

subplot(3,1,2);
plot(t_h, h_num, 'b-', 'LineWidth',2.0,'DisplayName','Numerical ||h||');
hold on;
plot(t_h, h_ref, 'r--', 'LineWidth',1.2,'DisplayName','Analytic Ref.');
xline(t_wait_h, 'k:', 'LineWidth',1.0, 'HandleVisibility','off');
grid on;
xlabel('Time [h]','FontWeight','bold'); ylabel('||h||  [km^2/s]','FontWeight','bold');
title('Specific Angular Momentum  ||h||');
legend('Location','best');

subplot(3,1,3);
plot(t_h, ecc_num, 'b-', 'LineWidth',2.0,'DisplayName','Numerical ||e||');
hold on;
plot(t_h, ecc_ref, 'r--', 'LineWidth',1.2,'DisplayName','Analytic Ref.');
xline(t_wait_h, 'k:', 'LineWidth',1.0, 'HandleVisibility','off');
grid on;
xlabel('Time [h]','FontWeight','bold'); ylabel('||e||  [-]','FontWeight','bold');
title('Eccentricity Magnitude  ||e||');
legend('Location','best');
sgtitle('Point 1.5 – LM Kepler Constants (Waiting + Transfer Phase)','FontSize',13);


%% =========================================================
%  Helper function
% ==========================================================
function s = verdict(err, tol)
    if err < tol
        s = 'PASS';
    else
        s = 'FAIL';
    end
end