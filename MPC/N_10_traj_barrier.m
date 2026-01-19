clear;
clc;


import casadi.*
starting_point = [1; 1; 0];

rob_radius = 0.2;

T = 0.1; %[s]
N = 10; % prediction horizon

v_max = 0.5; v_min = -v_max;
omega_max = pi/2; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)

[solver, args] = mpc_setup_barrier(N, 0.3);

obstacle1.center = [5 4];
obstacle1.radius = 0.5;
obstacle_sq_1_list = [2; 3; 3; 4];
obstacle_sq_2_list = [5.3; 5.5; 5.8; 6.5];

% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = starting_point;    % initial condition.
% xs = [1.5 ; 1.5 ; 0.0]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,2);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 30; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];
xs = [7; 7; 0];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
x_ref_jac = [];
y_ref_jac = [];
predicted_traj={};
current_time = 0;
ttt = 0;
last_u = zeros(N,2);
N_list = [];
while(norm((x0-xs),2) > 1e-1 && current_time < sim_tim) % new - condition for ending the loop

    % disp(x0)
    current_time = current_time + ttt;  %new - get the current time
    % args.p   = [x0;xs]; % set the values of the parameters vector
    %----------------------------------------------------------------------
    args.p(1:3) = x0; % initial condition of the robot posture
    for k = 1:N %new - set the reference to track
        x_ref = 7; y_ref = 7; theta_ref = 0;
        u_ref = 0; omega_ref = 0;
        args.p(5*k-1:5*k+1) = [x_ref, y_ref, theta_ref];
        args.p(5*k+2:5*k+3) = [u_ref, omega_ref];
        args.p(5*(N+1)-1:5*(N+1)+1) = [x_ref, y_ref, theta_ref];
    end
    x_ref_jac = [x_ref_jac x_ref];
    y_ref_jac = [y_ref_jac y_ref];
    
    %----------------------------------------------------------------------    
    % initial value of the optimization variables
    main_loop = tic;
    args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);

    status = solver.stats.return_status;  % 获取求解状态
    
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution

    if u(1, 1) < 1e-9 && u(1, 2) < 1e-9 && norm((x0-xs),2) > 0.2
        fprintf('陷入局部最优，尝试换初始值\n', status);
        u0 = zeros(N,2);        % two control inputs for each robot
        X0 = repmat(xs,1,N+1)'; % initialization of the states decision variables

        args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution
    end

    % 检查状态
    if ~strcmp(status, 'Solve_Succeeded')
        fprintf('求解失败，状态：%s 尝试换初始值\n', status);
        u0 = zeros(N,2);        % two control inputs for each robot
        X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

        args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
        status = solver.stats.return_status;  % 获取求解状态
        if ~strcmp(status, 'Solve_Succeeded')
            fprintf('求解失败，状态：%s 没救了，用上次输入\n', status);
        end

    end
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution

    predicted_traj{mpciter+1} = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    ttt = toc(main_loop);
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;

    if ttt > T/2
        disp(ttt);
    end
    last_x0 = x0;
    % disp(ttt);
    % Apply the control and shift the solution
    [t0, x0, u0] = shift_Runge_Kutta(ttt, t0, x0, last_u,f);
    % u0 = [u(1, :);u0(1:size(u,1)-1,:)];
    u0=u;
    last_u = u;
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    % % Shift trajectory to initialize the next step
    % disp([X0(1, 1:2)> x0(1:2)'])
    X0 = [x0';X0(2:end,:)];
    % X0 = repmat(xs,1,N+1)';
    % mpciter
    mpciter = mpciter + 1;

    obs1_size = sqrt(2);
    obs2_size = 1;
    obs3_size = 1.118;

    kk1 = 0.9;
    kk2 = 0.9*obs2_size/obs1_size;
    kk3 = 0.9*obs3_size/obs1_size;

    last_dist_1 = pointToRectangleDistance(last_x0(1), last_x0(2), obstacle_sq_1_list(1),obstacle_sq_1_list(2), obstacle_sq_1_list(3), obstacle_sq_1_list(4));
    curr_dist_1 = pointToRectangleDistance(x0(1), x0(2), obstacle_sq_1_list(1),obstacle_sq_1_list(2), obstacle_sq_1_list(3), obstacle_sq_1_list(4));

    last_dist_2 = sqrt((last_x0(1)-obstacle1.center(1))^2+(last_x0(2)-obstacle1.center(2))^2) - obstacle1.radius;
    curr_dist_2 = sqrt((x0(1)-obstacle1.center(1))^2+(x0(2)-obstacle1.center(2))^2) - obstacle1.radius;

    last_dist_3 = pointToRectangleDistance(last_x0(1), last_x0(2), obstacle_sq_2_list(1),obstacle_sq_2_list(2), obstacle_sq_2_list(3), obstacle_sq_2_list(4));
    curr_dist_3 = pointToRectangleDistance(x0(1), x0(2), obstacle_sq_2_list(1),obstacle_sq_2_list(2), obstacle_sq_2_list(3), obstacle_sq_2_list(4));
 
end



output1 = mpciter;
output2 = current_time;
output3 = mpciter / current_time;

xx_10_barrier = xx;

save("data_traj.mat", "xx_10_barrier","-append");


obstacle1.center = [5 4];
obstacle1.radius = 0.5;

obstacle_sq_1.duration  = [2 3];
obstacle_sq_1.width     = [4 3];

obstacle_sq_2.duration  = [5.3 5.8];
obstacle_sq_2.width     = [6.5 5.5];

figure(10);
hold on;
plot(xx(1,:), xx(2,:), "DisplayName", "Actual trajectory");
rectangle('Position',[obstacle1.center(1)-obstacle1.radius, obstacle1.center(2)-obstacle1.radius, 2*obstacle1.radius, 2*obstacle1.radius], 'Curvature',[1,1]);
rectangle('Position', [min(obstacle_sq_1.duration) min(obstacle_sq_1.width) abs(obstacle_sq_1.duration(2)-obstacle_sq_1.duration(1)) abs(obstacle_sq_1.width(2)-obstacle_sq_1.width(1))], 'Curvature', 0.3)
rectangle('Position', [min(obstacle_sq_2.duration) min(obstacle_sq_2.width) abs(obstacle_sq_2.duration(2)-obstacle_sq_2.duration(1)) abs(obstacle_sq_2.width(2)-obstacle_sq_2.width(1))], 'Curvature', 0.3)

xlim([0, 8])
ylim([0, 8])
title("mpciter = "+ mpciter + ", time = " + current_time)
axis equal
grid on;
legend






% figure('color','w')
% fig = figure(1);
% plot(t, N_list)
% xlabel('time(s)', Fontname='Times New Roman');
% ylabel('predictive horizon', Fontname='Times New Roman');
% export_fig 'predictive_horizon_fig' -eps
