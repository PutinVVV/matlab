close all; clc;

% System Parameters
N = 4;         % Number
l = 0.5;        % length
m = 1.0;        % Mass
J = 2/3*m*l^2;  % Moment of inert
ct = 0.1;       % Tangentialcoefficient
cn = 0.05;      % Normal coefficient

% Time parameters
t_start = 0;
t_end = 10;
dt = 0.01;
t = t_start:dt:t_end;

% Initial conditions
theta = zeros(N, 1);       % zero
theta_dot = zeros(N, 1);   % zero
x_cm = 0;                  % x cm
y_cm = 0;                  % y position cm
v_x = 0;                   % x speed
v_y = 0;                   % y speed

% Define A and D
A = zeros(N-1, N);
D = zeros(N-1, N);
for i = 1:N-1
    D(i,i) = 1;
    D(i,i+1) = -1;
    A(i,i) = 1;
    A(i,i+1) = 1;
end

% Create a structure
params.N = N;
params.l = l;
params.m = m;
params.J = J;
params.ct = ct;
params.cn = cn;
params.A = A;
params.D = D;

% Solve ODE45
y0 = [theta; theta_dot; x_cm; y_cm; v_x; v_y];
[t_sol, y_sol] = ode45(@(t,y) system_dynamics_with_cm(t,y,params), [t_start t_end], y0);

% Extract sol
theta_sol = y_sol(:, 1:N);
theta_dot_sol = y_sol(:, N+1:2*N);
x_cm_sol = y_sol(:, 2*N+1);
y_cm_sol = y_sol(:, 2*N+2);
v_x_sol = y_sol(:, 2*N+3);
v_y_sol = y_sol(:, 2*N+4);

figure;
subplot(2,1,1);
plot(t_sol, theta_sol);
title('Joint Angles');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend(arrayfun(@(i) sprintf('Link %d', i), 1:N, 'UniformOutput', false));
grid on;

subplot(2,1,2);
plot(t_sol, theta_dot_sol);
title('Angular Velocities');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend(arrayfun(@(i) sprintf('Link %d', i), 1:N, 'UniformOutput', false));
grid on;

% Animation
figure;
for i = 1:length(t_sol)
    x = zeros(N+1, 1);
    y = zeros(N+1, 1);
    
    % Add center of mass position
    x(1) = x_cm_sol(i);
    y(1) = y_cm_sol(i);
    
    for j = 1:N
        x(j+1) = x(j) + 2*l*cos(theta_sol(i,j));
        y(j+1) = y(j) + 2*l*sin(theta_sol(i,j));
    end
    
    clf;
    plot(x, y, 'b-o', 'LineWidth', 2);
    hold on;
    plot(x_cm_sol(i), y_cm_sol(i), 'r*', 'MarkerSize', 10);
    axis equal;
    xlim([min(x_cm_sol)-3*l*N max(x_cm_sol)+3*l*N]);
    ylim([min(y_cm_sol)-3*l*N max(y_cm_sol)+3*l*N]);
    grid on;
    title(sprintf('Time: %.2f s', t_sol(i)));
    drawnow;
    pause(0.01);
end

function [M_theta, W, R] = compute_matrices(theta, theta_dot, params)
    N = params.N;
    l = params.l;
    m = params.m;
    J = params.J;
    ct = params.ct;
    cn = params.cn;
    A = params.A;
    D = params.D;
    
    % Compute sin and cos matrices
    S_theta = diag(sin(theta));
    C_theta = diag(cos(theta));
    
    % Compute V and K 
    DDT = D * D';
    if det(DDT) == 0
        error('Matrix DDT is singular, cannot invert.');
    end
    DDT_inv = inv(DDT);
    V = A' * DDT_inv * A;
    K = A' * DDT_inv * D;
    
    % Compute M_theta
    M_theta = J*eye(N) + m*l^2*(S_theta*V*S_theta + C_theta*V*C_theta);
    
    % Compute W
    W = m*l^2*(S_theta*V*C_theta - C_theta*V*S_theta);
    
    % Compute friction
    Ct = ct*eye(N);
    Cn = cn*eye(N);
    
    R_x = Ct.*(C_theta.^2) + Cn.*(S_theta.^2);
    R_y = (Ct-Cn).*(S_theta.*C_theta);
    
    R = [R_x, R_y];
end

function dy = system_dynamics_with_cm(t, y, params)
    N = params.N;
    l = params.l;
    D = params.D;
    
    theta = y(1:N);
    theta_dot = y(N+1:2*N);
    x_cm = y(2*N+1);
    y_cm = y(2*N+2);
    v_x = y(2*N+3);
    v_y = y(2*N+4);
    
    % Compute matrices
    [M_theta, W, R] = compute_matrices(theta, theta_dot, params);
    
    % Compute K matrix
    A = params.A;
    DDT = D * D';
    if det(DDT) == 0
        error('Matrix DDT is singular, cannot invert.');
    end
    DDT_inv = inv(DDT);
    K = A' * DDT_inv * D;
    
    % Extract friction forces
    fR_x = R(:, 1);
    fR_y = R(:, 2);
    
    % Создаем управляющий сигнал для каждого звена
    u = zeros(N, 1);
    for i = 1:N
        u(i) = 0;
    end
    
    % Compute sin and cos matrices
    S_theta = diag(sin(theta));
    C_theta = diag(cos(theta));
    
    % Вычисляем правую часть уравнения
    RHS = D'*u(1:N-1);
    
    % Вычисляем левую часть уравнения
    term1 = W * (theta_dot.^2);
    term2 = l * (S_theta * K * fR_x - C_theta * K * fR_y);
    LHS = term1 - term2;
    
    % Решаем уравнение M_theta*theta_two_dot = RHS - LHS
    if det(M_theta) == 0
        error('Matrix M_theta is singular, cannot solve for theta_ddot.');
    end
    theta_ddot = M_theta \ (RHS - LHS);
    
    % Вычисляем ускорение центра масс
    a_x = sum(fR_x) / (N * params.m);
    a_y = sum(fR_y) / (N * params.m);
    
    dy = [theta_dot; theta_ddot; v_x; v_y; a_x; a_y];
end