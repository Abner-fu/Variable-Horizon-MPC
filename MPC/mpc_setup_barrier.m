function [solver, args] = mpc_setup_barrier(N, gamma)

import casadi.*

T = 0.1; %[s]
rob_radius = 0.2;

v_max = 0.5; v_min = -v_max;
omega_max = pi/2; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
%P = SX.sym('P',n_states + n_states);
P = SX.sym('P',n_states*2 + N*(n_states+n_controls));
% parameters (which include the initial state and the reference along the
% predicted trajectory (reference states and reference controls))

X = SX.sym('X',n_states,(N+1));
% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 10;Q(2,2) = 10;Q(3,3) = 0.1; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)
QN = zeros(3,3); QN(1,1) = 1000;QN(2,2) = 1000;QN(3,3) = 1;

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
for k = 1:N
    st = X(:,k);  con = U(:,k);
    %obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
    obj = obj+(st-P(5*k-1:5*k+1))'*Q*(st-P(5*k-1:5*k+1)) + ...
              (con-P(5*k+2:5*k+3))'*R*(con-P(5*k+2:5*k+3)) ; % calculate obj
    % the number is (n_states+n_controls)
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + T/2*k1, con); % new
    k3 = f(st + T/2*k2, con); % new
    k4 = f(st + T*k3, con); % new
    st_next_RK4=st +T/6*(k1 +2*k2 +2*k3 +k4); % new    
    % f_value = f(st,con);
    % st_next_euler = st+ (h*f_value);
    % g = [g;st_next-st_next_euler]; % compute constraints
    g = [g;st_next-st_next_RK4]; % compute constraints % new
end
obj = obj + (X(:,N+1)-P(5*(N+1)-1:5*(N+1)+1))'*QN*(X(:,N+1)-P(5*(N+1)-1:5*(N+1)+1));

% Add constraints for collision avoidance
obstacle1.center = [5 4];
obstacle1.radius = 0.5;
obstacle2.center = [4 5];
obstacle2.radius = 1;
obstacle3.center = [6 2];
obstacle3.radius = 0.3;


% for k = 1:N+1   % obstacle1 constraints
%     g = [g ; -((X(1,k)-obstacle1.center(1))^2+(X(2,k)-obstacle1.center(2))^2) + (obstacle1.radius + rob_radius)^2];
% end
% for k = 1:N+1   % obstacle1 constraints
%     g = [g ; -((X(1,k)-obstacle2.center(1))^2+(X(2,k)-obstacle2.center(2))^2) + (obstacle2.radius + rob_radius)^2];
% end
% for k = 1:N+1   % obstacle1 constraints
%     g = [g ; -((X(1,k)-obstacle3.center(1))^2+(X(2,k)-obstacle3.center(2))^2) + (obstacle3.radius + rob_radius)^2];
% end

for k = 1:N   % cbf constraints
    hh1 = ((X(1,k)-obstacle1.center(1))^2+(X(2,k)-obstacle1.center(2))^2) - ((obstacle1.radius + rob_radius)^2);
    hh2 = ((X(1,k+1)-obstacle1.center(1))^2+(X(2,k+1)-obstacle1.center(2))^2) - ((obstacle1.radius + rob_radius)^2);
    g = [g ; (1-gamma)*hh1-hh2];
end

% for k = 1:N   % cbf constraints
%     h3 = ((X(1,k)-obstacle2.center(1))^2+(X(2,k)-obstacle2.center(2))^2) - ((obstacle2.radius + rob_radius)^2);
%     h4 = ((X(1,k+1)-obstacle2.center(1))^2+(X(2,k+1)-obstacle2.center(2))^2) - ((obstacle2.radius + rob_radius)^2);
%     g = [g ; (1-gamma)*h3-h4];
% end

% for k = 1:N   % cbf constraints
%     h5 = ((X(1,k)-obstacle3.center(1))^2+(X(2,k)-obstacle3.center(2))^2) - ((obstacle3.radius + rob_radius)^2);
%     h6 = ((X(1,k+1)-obstacle3.center(1))^2+(X(2,k+1)-obstacle3.center(2))^2) - ((obstacle3.radius + rob_radius)^2);
%     g = [g ; (1-gamma)*h5-h6];
% end
obstacle_sq_1.duration  = [2 3];
obstacle_sq_1.width     = [4 3];

obstacle_sq_2.duration  = [5.3 5.8];
obstacle_sq_2.width     = [6.5 5.5];
% for k = 1:N+1   % obstacle1 constraints
%     g = [g ; min([X(1,k)-(obstacle_sq_1.duration(1)-rob_radius), ...
%         (obstacle_sq_1.duration(2)+rob_radius)-X(1,k),...
%         X(2,k)-(obstacle_sq_1.width(2)-rob_radius), ...
%         (obstacle_sq_1.width(1)+rob_radius)-X(2,k)])];
% end
for k = 1:N   % obstacle1 constraints
    hs1 = min([X(1,k)-(obstacle_sq_1.duration(1)-rob_radius), ...
        (obstacle_sq_1.duration(2)+rob_radius)-X(1,k),...
        X(2,k)-(obstacle_sq_1.width(2)-rob_radius), ...
        (obstacle_sq_1.width(1)+rob_radius)-X(2,k)]);
    hs2 = min([X(1,k+1)-(obstacle_sq_1.duration(1)-rob_radius), ...
        (obstacle_sq_1.duration(2)+rob_radius)-X(1,k+1),...
        X(2,k+1)-(obstacle_sq_1.width(2)-rob_radius), ...
        (obstacle_sq_1.width(1)+rob_radius)-X(2,k+1)]);
    g = [g ; hs2-(1-gamma)*hs1];
end

for k = 1:N   % obstacle1 constraints
    hs3 = min([X(1,k)-(obstacle_sq_2.duration(1)-rob_radius), ...
        (obstacle_sq_2.duration(2)+rob_radius)-X(1,k),...
        X(2,k)-(obstacle_sq_2.width(2)-rob_radius), ...
        (obstacle_sq_2.width(1)+rob_radius)-X(2,k)]);
    hs4 = min([X(1,k+1)-(obstacle_sq_2.duration(1)-rob_radius), ...
        (obstacle_sq_2.duration(2)+rob_radius)-X(1,k+1),...
        X(2,k+1)-(obstacle_sq_2.width(2)-rob_radius), ...
        (obstacle_sq_2.width(1)+rob_radius)-X(2,k+1)]);
    g = [g ; hs4-(1-gamma)*hs3];
end

% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 20000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-4;
opts.ipopt.mu_init = 1e-4;
opts.ipopt.acceptable_obj_change_tol = 1e-4;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg(1:3*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:3*(N+1)) = 0;  % 1e-20   % Equality constraints

args.lbg(3*(N+1)+1 : 3*(N+1)+ 3*N) = -inf; % inequality constraints
args.ubg(3*(N+1)+1 : 3*(N+1)+ 3*N) = 0; % inequality constraints

args.lbx(1:3:3*(N+1),1) = -1; %state x lower bound % new - adapt the bound
args.ubx(1:3:3*(N+1),1) = 8; %state x upper bound  % new - adapt the bound
args.lbx(2:3:3*(N+1),1) = -1; %state y lower bound
args.ubx(2:3:3*(N+1),1) = 8; %state y upper bound
args.lbx(3:3:3*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:3:3*(N+1),1) = inf; %state theta upper bound

args.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_min; %v lower bound
args.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_max; %v upper bound
args.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_min; %omega lower bound
args.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_max; %omega upper bound


end

