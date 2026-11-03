clear all;
close all;
clc;

% Part 4 - Kinematics and force propagation of the 5R berthing manipulator.
%
% The arm is modelled as an open-chain multibody system with 5 revolute joints.
% Forward kinematics is used to compute the pose of each link and of the
% end-effector. Differential kinematics is then formulated through the
% System Jacobian matrix N, allowing the propagation of velocities,
% wrenches and generalized forces throughout the structure.

%% ============================================
% PART 4.2 - FORWARD KINEMATICS
% ============================================
% A coordinate frame is attached to each link.
% The pose of link Li with respect to the inertial frame is obtained
% through successive Homogeneous Transformation Matrices (HTMs).

L = [2.5, 2, 1.5, 1.2, 0.8];    % Arm Lenght [m]

% Four representative configurations are analysed:
% 1) Trivial configuration
% 2) Non-trivial configuration
% 3) Minimum torque configuration
% 4) Maximum torque configuration

theta = [0, 0, 0, 0, 0;         % Trivial Configuration
         0, pi/4, 0, pi, 0;     % Non-Trival Configuration
         0, pi/2, 0, pi/2, 0;   % Configuration expected to minimise joint torques
         0, 0, 0, pi/2, 0];     % Configuration expected to maximise joint torques
         
syms theta1 theta2 theta3 theta4 theta5 real
syms L1 L2 L3 L4 L5 real

T_EE_sym= HTM([theta1;theta2;theta3;theta4;theta5],[L1,L2,L3,L4,L5]);
T_EE_2D = subs(T_EE_sym, [theta1, theta3, theta5], [0, 0, 0]);
T_EE_2D = simplify(T_EE_2D);

disp('Homogeneous Transformation Matrix (HTM) - 2D');
disp(T_EE_2D);

for f = 1:4
    theta_current = theta(f,:);

    [T_EE_num,T01,T02,T03,T04,T05,TEE] = HTM(theta_current,L);

    P1 = T01(1:3, 4);         % J1
    P2 = T02(1:3, 4);         % J2
    P3 = T03(1:3, 4);         % J3
    P4 = T04(1:3, 4);         % J4
    P5 = T05(1:3, 4);         % J5
    PEE = TEE(1:3, 4);        % End - Effector

    X = [P1(1), P2(1), P3(1), P4(1), P5(1), PEE(1)];
    Y = [P1(2), P2(2), P3(2), P4(2), P5(2), PEE(2)];
    Z = [P1(3), P2(3), P3(3), P4(3), P5(3), PEE(3)];

    figure;
    plot3(X, Y, Z, '-o', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'Color', 'b');
    hold on;
    grid on;
    axis equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');

    if f==1
    title('Trival Configuration - \theta_{2} = 0 and \theta_{4} = 0');
    elseif f==2
    title('Non-Trival Configuration - \theta_{2} = \pi/4 and \theta_{4} = \pi');
    elseif f==3
    title('Configuration that minimises the torque - \theta_{2} = \pi/2 and \theta_{4} = \pi/2');
    else 
    title('Configuration that maximises the torque - \theta_{2} = \0 and \theta_{4} = \pi/2');

    end

plot3(PEE(1), PEE(2), PEE(3), '*', 'MarkerSize', 6, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y');
patch([-3 3 3 -3], [-3 -3 3 3], [0 0 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
view(3);

end 


% HTM Function
function [T_EE,T01,T02,T03,T04,T05,TEE] = HTM(theta,L)

T_0_1 = [cos(theta(1)), -sin(theta(1)), 0, 0;
         sin(theta(1)), cos(theta(1)), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];

T_1_2 = [1, 0, 0, 0;
    0, 0, -1, 0;
    0, 1, 0, L(1);
    0, 0, 0, 1] * ...
    [cos(theta(2)), -sin(theta(2)), 0, 0;
    sin(theta(2)), cos(theta(2)), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1] * eye(4);

T_2_3 = [1, 0, 0, L(2);
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1] * ...
    [cos(theta(3)), -sin(theta(3)), 0, 0;
    sin(theta(3)), cos(theta(3)), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1] * eye(4);

T_3_4 = [1, 0, 0, L(3);
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1] * ...
    [cos(theta(4)), -sin(theta(4)), 0, 0;
    sin(theta(4)), cos(theta(4)), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1] * eye(4);

T_4_J5 = [1, 0, 0, 0;
    0, 0, -1, -L(4);
    0, 1, 0, 0;
    0, 0, 0, 1] * ...
    [cos(theta(5)), -sin(theta(5)), 0, 0;
    sin(theta(5)), cos(theta(5)), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1];
T_J5_EE = [1, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, 1, L(5);
    0, 0, 0, 1];

% Compute the overall HTM

T_4_5 = T_4_J5*T_J5_EE;
T_EE = T_0_1*T_1_2*T_2_3*T_3_4*T_4_5;

T01 = T_0_1;
T02 = T01*T_1_2;
T03 = T02*T_2_3;
T04 = T03*T_3_4;
T05 = T04*T_4_J5;
TEE = T05*T_J5_EE;


end

%% ============================================
% PART 4.3 - DIFFERENTIAL KINEMATICS
% ============================================

% Twist propagation is described through:
% p_i  : twist propagation vectors
% B_ij : twist propagation matrices
% The system Jacobian is assembled as: N = N_l N_d

[N,Nl,Nd] = JacobianSystem(theta_current,L);

[rows, columns] = size(N);   
fprintf('Extracted dimensions of N: %d x %d\n',rows, columns);

function [N,Nl,Nd] = JacobianSystem(theta,L)

n = 5;

[~,T01,T02,T03,T04,T05,TEE] = HTM(theta,L);

% Rotation matrices (3x3)
R01 = T01(1:3,1:3);
R02 = T02(1:3,1:3);
R03 = T03(1:3,1:3);
R04 = T04(1:3,1:3);
R05 = T05(1:3,1:3);

% ============================================
% Link CoM positions
% ============================================

% Local CoM vectors

g1_local = [0; 0; L(1)/2];
g2_local = [L(2)/2; 0; 0];
g3_local = [L(3)/2; 0; 0];
g4_local = [0; -L(4)/2; 0];
g5_local = [0; 0; L(5)/2];

% CoM vectors expressed in inertial frame

g1 = R01 * g1_local;
g2 = R02 * g2_local;
g3 = R03 * g3_local;
g4 = R04 * g4_local;
g5 = R05 * g5_local;

% Joint positions expressed in inertial frame (3x1)
P0 = [0;0;0];
P1 = T01(1:3,4);
P2 = T02(1:3,4);
P3 = T03(1:3,4);
P4 = T04(1:3,4);
P5 = T05(1:3,4);

% Absolute CoM positions

C0 = [0;0;0];
C1 = P1 + g1;
C2 = P2 + g2;
C3 = P3 + g3;
C4 = P4 + g4;
C5 = P5 + g5;

% Joint axes:
% the third column of each rotation matrix corresponds
% to the local z-axis expressed in the inertial frame.
% Each z_i is a (3x1) vector.

z0 = [0;0;1];
z1 = R01(:,3);
z2 = R02(:,3);
z3 = R03(:,3);
z4 = R04(:,3);
z5 = R05(:,3);

% ============================================
% Twist Propagation Vectors p_i
% ============================================

% Each p_i is a (6x1) vector:
%
% p_i = [ z_i ;
%         z_i x g_i ]
%
% upper part  -> angular velocity contribution
% lower part  -> linear velocity contribution

p1 = [z1;
      cross(z1,g1)];

p2 = [z2;
      cross(z2,g2)];

p3 = [z3;
      cross(z3,g3)];

p4 = [z4;
      cross(z4,g4)];

p5 = [z5;
      cross(z5,g5)];

% ============================================
% Skew matrix function
% ============================================

skew = @(v)[  0   -v(3)  v(2);
             v(3)   0   -v(1);
            -v(2) v(1)    0 ];

% ============================================
% Twist Propagation Matrices B_ij
% ============================================

% Each B_ij is a (6x6) matrix:
%
% B_ij = [ I3      0
%          skew(r) I3 ]
%
% where r = C_i - C_j

B10 = [eye(3) zeros(3);
       skew(C1-C0) eye(3)];

B21 = [eye(3) zeros(3);
       skew(C2-C1) eye(3)];

B32 = [eye(3) zeros(3);
       skew(C3-C2) eye(3)];

B43 = [eye(3) zeros(3);
       skew(C4-C3) eye(3)];

B54 = [eye(3) zeros(3);
       skew(C5-C4) eye(3)];

% ============================================
% LOWER TRIANGULAR MATRIX N_l
% ============================================

% N_l dimensions: 6(n+1) x 6(n+1)
% For n = 5: N_l ∈ R^(36x36)

Nl = zeros(6*(n+1), 6*(n+1),'like',theta);

% Identity blocks on diagonal
%
% Each body contributes with a 6x6 block.
%
% rows = 6*(i-1)+1 : 6*i
% Example:
% i=1 -> rows 1:6
% i=2 -> rows 7:12
% i=3 -> rows 13:18
% ...

for i = 1:n+1
    
    rows = 6*(i-1)+1 : 6*i;
    
    Nl(rows,rows) = eye(6,'like',theta);
    
end

% Store all B matrices inside a 3D array
%
% B(:,:,1) = B10
% B(:,:,2) = B21
% ...

B(:,:,1) = B10;
B(:,:,2) = B21;
B(:,:,3) = B32;
B(:,:,4) = B43;
B(:,:,5) = B54;

% Lower triangular propagation
%
% The element (i,j) of N_l contains the propagation
% from body j to body i.
%
% Example:
%
% N_l(13:18,1:6) = B32*B21*B10
%
% which propagates the twist from base to link 3.

for i = 2:n+1
    
    for j = 1:i-1
        
        Bij_product = eye(6,'like',theta);
        
        for k = j:i-1
            
            Bij_product = B(:,:,k) * Bij_product;
            
        end
        
        rows = 6*(i-1)+1 : 6*i;
        cols = 6*(j-1)+1 : 6*j;
        
        Nl(rows,cols) = Bij_product;
        
    end
    
end

% ============================================
% BLOCK DIAGONAL MATRIX N_d
% ============================================

% N_d dimensions: 6(n+1) x (6+n)
% For n = 5: N_d ∈ R^(36x11)

Nd = zeros(6*(n+1),6+n,'like', theta);

% Base block:
% The first 6 columns correspond to the base twist.

Nd(1:6,1:6) = eye(6,'like',theta);

% List of all p_i vectors
% P_twist dimensions: 6 x n

P_twist = [p1 p2 p3 p4 p5];

% Joint blocks
% Each joint contributes with ONE column because
% every revolute joint has 1 DoF.
%
% Example:
% Joint 1 -> column 7
% Joint 2 -> column 8
% ...

for i = 1:n
    
    % row_start / row_end select the rows associated to body i.
    row_start = 6*i + 1;
    row_end   = 6*(i+1);
    
    col = 6 + i;
    
    Nd(row_start:row_end,col) = P_twist(:,i);
    
end

% ============================================
% SYSTEM JACOBIAN
% ============================================

% Final System Jacobian: N = N_l * N_d
% Dimensions: (36x36)*(36x11) = (36x11)

N = Nl*Nd;
end

%%
% ============================================
% PART 4.4 - CONTACT FORCE PROPAGATION
% ============================================

% CONTACT WRENCH ON LINK L5
%
% We define a contact wrench applied ONLY at the end-effector.
%
% The chosen wrench is a lateral contact force along inertial x-axis with
% no external moment.
% Physically this models the contact force generated during docking, acting
% along the final approach direction.

w5 = [0;     % nx [Nm]
      0;     % ny [Nm]
      0;     % nz [Nm]
      0;     % fx [N]
      0;     % fy [N]
      100];    % fz [N]

% ============================================
% FULL SYSTEM WRENCH VECTOR
% ============================================

% Total dimensions: w ∈ R^(6*(n+1))
%
% Only the last 6 elements are non-zero because
% the wrench is applied only on the final link.
n = 5;
w = zeros(6*(n+1),1);

% Insert wrench on link L5
w(end-5:end) = w5;

% ============================================
% GENERALISED FORCE VECTOR
% ============================================

% System equation:
tau = N' * w;

% tau = [tau0 ; taum]
% tau0 ∈ R^6  -> base reactions
% taum ∈ R^n -> joint torques

tau0 = tau(1:6);
taum = tau(7:end);

% ============================================
% RESULTS
% ============================================

% ============================================
% RESULTS FOR ALL CONFIGURATIONS
% ============================================

for f = 1:size(theta,1)

    theta_current = theta(f,:);

    % System Jacobian for current configuration
    [N,~,~] = JacobianSystem(theta_current,L);

    % Generalised forces
    tau = N' * w;

    tau0 = tau(1:6);
    taum = tau(7:end);

    fprintf('\n========================================\n')

    switch f
        case 1
            fprintf('TRIVIAL CONFIGURATION\n')
        case 2
            fprintf('NON-TRIVIAL CONFIGURATION\n')
        case 3
            fprintf('MINIMUM TORQUE CONFIGURATION\n')
        case 4
            fprintf('MAXIMUM TORQUE CONFIGURATION\n')
    end

    fprintf('========================================\n')

    fprintf('\nBase reaction wrench tau0 = [n0 ; f0]\n')
    fprintf('nx = %8.2f Nm\n',tau0(1))
    fprintf('ny = %8.2f Nm\n',tau0(2))
    fprintf('nz = %8.2f Nm\n',tau0(3))
    fprintf('fx = %8.2f N\n',tau0(4))
    fprintf('fy = %8.2f N\n',tau0(5))
    fprintf('fz = %8.2f N\n',tau0(6))
    
    fprintf('\nJoint torques taum\n')
    for i = 1:n
        fprintf('tau%d = %8.2f Nm\n',i,taum(i)) 
    end

end


%%
% ============================================
% BONUS POINT - INVERSE KINEMATICS
% ============================================

% Initial joint configuration 
theta_start = [0; 0; 0; pi/2; 0];
theta_k = theta_start;

% Target end-effector position 
x_target_pos = [4; 0; 1]; %[m]


% Algorithm parameters
tol = 1e-4;       % Convergence tolerance on position error norm
lambda = 0.05;    % Damping factor for Damped Least Squares (DLS)
K = 0.5;          % Step gain
max_iter = 10000; % Maximum number of interations

% Initialisation
iter = 0;
error_norm = 1;
error_history = [];

% Iterative IK loop - Position only
while (error_norm>tol) && (iter < max_iter)

    % Forward Kinematics: Compute current end-effector pose
    T_EE_current = HTM(theta_k,L);
    pose = T_EE_current(1:3,4);
    
    % Position error in Cartisian space
    delta_pos = x_target_pos - pose;

    % Update error norm and store history for convergence plot
    error_norm = norm(delta_pos);
    error_history(end+1) = error_norm;

    % Compute geomtric Jacobian
    J_full = GeometricJacobian(theta_k,L);
    J = J_full(1:3,:); % Extract position-only part (3x5)

    % Damped Least Sqaures pseudo-inverse: J* = J^T (J J^T + lambda^2 I)^{-1}
    J_DLS = J' / (J*J' + lambda^2*eye(3));
 
    % Joint angle update
    theta_k = theta_k+K*J_DLS*delta_pos;

    iter = iter + 1;
end

% Store converged solution
theta_star = theta_k;

fprintf('\n========================================\n')
fprintf('INVERSE KINEMATICS\n')
fprintf('========================================\n')

disp('Converged Theta* (degrees):');
disp(rad2deg(wrapToPi(theta_star)));
disp(['Final Error Norm: ', num2str(error_norm)]);

% Forward kinematics at the converged configuration
% Extract joint and end-effector positions for plotting
[T_EE_num,T01,T02,T03,T04,T05,TEE] = HTM(theta_star,L);

P1 = T01(1:3, 4);         % J1
P2 = T02(1:3, 4);         % J2
P3 = T03(1:3, 4);         % J3
P4 = T04(1:3, 4);         % J4
P5 = T05(1:3, 4);         % J5
PEE = TEE(1:3, 4);        % End - Effector

% Assemble coordinate arrays
X = [P1(1), P2(1), P3(1), P4(1), P5(1), PEE(1)];
Y = [P1(2), P2(2), P3(2), P4(2), P5(2), PEE(2)];
Z = [P1(3), P2(3), P3(3), P4(3), P5(3), PEE(3)];

% Plot the manipulator at the IK solution configuration
figure;
plot3(X, Y, Z, '-o', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'Color', 'b');
hold on;
grid on;
axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('IK Solution - Converged Configuration \theta^*');

% Mark end-effector position
plot3(PEE(1), PEE(2), PEE(3), '*', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y');

% Mark target position
plot3(x_target_pos(1), x_target_pos(2), x_target_pos(3), '*g', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Reference ground plane
patch([-3 3 3 -3], [-3 -3 3 3], [0 0 0 0], [0.8 0.8 0.8], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');

legend('Manipulator links', 'End-effector', 'Target', 'Location', 'best');
view(3)


% GeometricJacobian - Computes the geometric Jacobian of the 5R manipulator
function J = GeometricJacobian(theta,L)

% Compute forward kinematics 
[~,T01,T02,T03,T04,T05,TEE] = HTM(theta,L);

% End-Effector Position in inertial frame
p_EE = TEE(1:3,4);

% Joint rotation axes expressed in the inertial frame
z0 = [0;0;1];    
z1 = T01(1:3, 3);
z2 = T02(1:3, 3);
z3 = T03(1:3, 3);
z4 = T04(1:3, 3);

% Joint origins expressed in the inertial frame
p0 = [0;0;0];
p1 = T01(1:3, 4);
p2 = T02(1:3, 4);
p3 = T03(1:3, 4);
p4 = T04(1:3, 4);

% Geometric Jacobian for Revolute Joint
J1 = [cross(z0,p_EE-p0); z0];
J2 = [cross(z1, p_EE-p1); z1];
J3 = [cross(z2, p_EE-p2); z2];
J4 = [cross(z3, p_EE-p3); z3];
J5 = [cross(z4, p_EE -p4); z4];

% Assemble the full 6x5 geometric jacobian
J = [J1,J2,J3,J4,J5];

end


%% Joint-Space Trajectory Planning via Cubic Polynomial Interpolation
T = 1800;  %[s] = 30 min
dt = 0.1;  % Time step [s]
t = 0:dt:T;
num_steps = length(t);

% Pre-allocate trajectory arrays
pos = zeros(5, num_steps);
vec = zeros(5, num_steps);
acc = zeros(5, num_steps);

% Compute cubic polynomial coefficents and evaluate trajectory for each
% joint
for i = 1:5

% Cubic polynomial coefficients    
a0 = theta_start(i);
a1 = 0;
a2 = 3*(theta_star(i) - theta_start(i))/ T^2;
a3 = -2*(theta_star(i)-theta_start(i)) / T^3;

% Joint angle trajectory 
pos(i,:) = a0 + a1*t + a2*t.^2 + a3*t.^3;

% Joint velocity (first derivative of cubic polynomial)
vec(i,:) = a1 + 2*a2*t + 3*a3*t.^2;

% Joint acceleration (second derivative - linear in time)
acc(i,:) = 2*a2 + 6*a3*t;

end

% Plot 1: Joint angle profiles
figure;
plot(t, rad2deg(pos), 'LineWidth', 2);
title('Position of the Joints');
xlabel('Time [s]');
ylabel('Angular Position [degrees]');
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', 'Location', 'best');
grid on;

% Plot 2: Joint velocity profiles
figure;
plot(t, rad2deg(vec), 'LineWidth', 2);
title('Joint Velocity Profile');
xlabel('Time [s]');
ylabel('Angular Velocity [deg/s]');
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', 'Location', 'best');
grid on;

% Plot 3: Joint acceleration profiles
figure;
plot(t, rad2deg(acc), 'LineWidth', 2);
title('Joints Accelerations Profile');
xlabel('Time [s]');
ylabel('Angular Acceleration [deg/s^2]');
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4', '\theta_5', 'Location', 'best');
grid on;

