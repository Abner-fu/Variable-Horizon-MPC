function [outputArg1,outputArg2] = fun_name(distance_obs, predicted_horizon_long)
predicted_horizon_long = round(predicted_horizon_long);
addpath('D:\casadi\casadi-3.6.6-windows64-matlab2018b')
import casadi.*

rob_radius = 0.3;

T = 0.1; %[s]
N = 2; % prediction horizon

v_max = 0.5; v_min = -v_max;
omega_max = pi/2; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)

[solver, args] = mpc_setup(N);



obstacle_sq_1.duration  = [3 4];
obstacle_sq_1.width     = [5 4];

obstacle_sq_1_list = [3; 4; 4; 5];

% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0 ; 3 ; 0.0];    % initial condition.
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
        % fprintf('陷入局部最优，尝试换初始值\n', status);
        u0 = zeros(N,2);        % two control inputs for each robot
        X0 = repmat(xs,1,N+1)'; % initialization of the states decision variables

        args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
        u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution
    end

    % 检查状态
    if ~strcmp(status, 'Solve_Succeeded')
        % fprintf('求解失败，状态：%s 尝试换初始值\n', status);
        u0 = zeros(N,2);        % two control inputs for each robot
        X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

        args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
        sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    
        status = solver.stats.return_status;  % 获取求解状态
        if ~strcmp(status, 'Solve_Succeeded')
            % fprintf('求解失败，状态：%s 没救了\n', status);
            outputArg1 = 10000;
            outputArg2 = 10000;
            return;
        end

    end
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution

    predicted_traj{mpciter+1} = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    ttt = toc(main_loop);
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;

   
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

    % last_dist = pointToRectangleDistance(last_x0(1), last_x0(2), obstacle_sq_1_list(1),obstacle_sq_1_list(2), obstacle_sq_1_list(3), obstacle_sq_1_list(4));
    curr_dist = pointToRectangleDistance(x0(1), x0(2), obstacle_sq_1_list(1),obstacle_sq_1_list(2), obstacle_sq_1_list(3), obstacle_sq_1_list(4));

    close_obs_flag = (curr_dist <= distance_obs && x0(1) < 3 && x0(2) < 5);
    if close_obs_flag && N ~= predicted_horizon_long
        N = predicted_horizon_long;
        [solver, args] = mpc_setup(N);
   
        X0 = [X0; repmat(X0(end,:), N+1-size(X0, 1), 1)];
        u0 = [u0; repmat(u0(end,:), N-size(u0, 1), 1)];
    end

    if ~close_obs_flag && N~=2
        N = 2;
        [solver, args] = mpc_setup(N);

        X0 = X0(1:N+1, :);
        u0 = u0(1:N, :);
    end

end

outputArg1 = 200 - (mpciter / current_time);
outputArg2 = compute_traj_dist(xx');

end

