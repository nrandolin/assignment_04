clear;

DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
DormandPrince.A = [0,0,0,0,0,0,0;
1/5, 0, 0, 0,0,0,0;...
3/40, 9/40, 0, 0, 0, 0,0;...
44/45, -56/15, 32/9, 0, 0, 0,0;...
19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

t = 2;
t_list = linspace(0,0.5,100);


t0 = 0;          % Start time
tf = 40;        % End time (replace with the desired final time)

tspan = [t0, tf];  % Full integration interval

X0 = 1;      % Initial conditions
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 6;
y0 = 10;
dxdt0 = 0;
dydt0 = 1.5;
my_rate_func = @(t_in,V_in) gravity_rate_func(t_in,V_in,orbit_params);

V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,40,500);
X_list_true = compute_planetary_motion(t_range,V0,orbit_params);

h_list = linspace(0.001, 0.12, 100);

%[XB, num_evals] = explicit_RK_step(@rate_func01,t,XA,h,BT_struct);
global_error_list = [];
h_avg_list = [];
num_evals_list = [];
for i = 1:length(h_list)
    h = h_list(i);
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
    (my_rate_func,tspan,V0,h,DormandPrince);
    X_numerical = X_list(end, :)';
    X_analytical = compute_planetary_motion(tspan(2),V0,orbit_params);
    global_error = norm(X_numerical - X_analytical);
    global_error_list_fixed = [global_error_list, global_error];
    h_avg_list_fixed = [h_avg_list, h_avg];
    num_evals_list_fixed = [num_evals_list, num_evals];
end

% plot planet trajectory
figure()
hold on
plot(X_list_true(:,1), X_list_true(:,2),'b','linewidth',2) 
plot(X_list(:,1), X_list(:,2), 'r--','linewidth',2)
title("Approximated vs. True Planet Trajectory")
legend("True Solution", "Approximated Solution")

% global error vs. avg step size
figure()
loglog(h_avg_list_fixed, global_error_list_fixed, 'b', 'LineWidth', 2)
title("Global Error vs. Average Step Size - Fixed Step")
xlabel("Average Step Size")
ylabel("Global Error")

% global error vs. num evals
figure()
loglog(num_evals_list_fixed, global_error_list_fixed, 'b', 'LineWidth', 2)
title("Global Error vs. Number of Evaluations - Fixed Step")
xlabel("Number of Evaluations")
ylabel("Global Error")

%% Heun's Method
heun_struct.A = [0, 0; 1, 0];
heun_struct.B = [0.5, 0.5];
heun_struct.C = [0; 1];
h = 0.001;
t0 = 0;          % Start time
tf = 40;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval  

[t_list_heun,X_list_heun,h_avg_heun, num_evals_heun] = explicit_RK_fixed_step_integration ...
(my_rate_func,tspan,V0,h,heun_struct);

figure()
hold on
plot(X_list_true(:,1), X_list_true(:,2),'b','linewidth',2) 
plot(X_list_heun(:,1), X_list_heun(:,2), 'r--','linewidth',2)
xlabel("X")
ylabel("Y")
title("Heun's Method Approximated Solution")
legend("Heun's Method", "True Solution")

% Local Error Calculations
h_list = logspace(-5,1,100); 

local_error_heun = [];
X_true = [];
difference = [];

t_ref = 4.49;
XA = compute_planetary_motion(t_ref, V0, orbit_params);

for i = 1:length(h_list)
   h_ref = h_list(i);
   [XB_heun,~] = explicit_RK_step(my_rate_func,t_ref,XA,h_ref,heun_struct);

   X_analytical = compute_planetary_motion(t_ref+h_ref, V0, orbit_params);
    
   X_true = [X_true, X_analytical];
   local_error_heun = [local_error_heun, norm(XB_heun - X_analytical)];
   difference = [difference, norm(X_analytical-compute_planetary_motion(t_ref, V0, orbit_params))];

end
figure();
loglog(h_list, local_error_heun, 'b', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error Heun's Method")
legend("Heun Method's", "Difference")

%% Ralston's Method - Third Order
ralston_struct.A = [0, 0, 0; 0.5, 0, 0; 0, 0.75, 0];
ralston_struct.B = [2/9, 1/3, 4/9];
ralston_struct.C = [0; 0/5; 0.75];
h = 0.001;
t0 = 0;          % Start time
tf = 40;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval  

[t_list_ralston,X_list_ralston,h_avg_ralston, num_evals_ralston] = explicit_RK_fixed_step_integration ...
(my_rate_func,tspan,V0,h,ralston_struct);

figure()
hold on
plot(X_list_true(:,1), X_list_true(:,2),'b','linewidth',2) 
plot(X_list_ralston(:,1), X_list_ralston(:,2), 'r--','linewidth',2)
xlabel("X")
ylabel("Y")
title("Ralston's Method Approximated Solution")
legend("Ralston's Method", "True Solution")

% Local Error Calculations
h_list = logspace(-5,1,100); 

local_error_ralston = [];
X_true = [];
difference = [];

t_ref = 4.49;
XA = compute_planetary_motion(t_ref, V0, orbit_params);

for i = 1:length(h_list)
   h_ref = h_list(i);
    
   [XB_ralston,~] = explicit_RK_step(my_rate_func,t_ref,V0,h_ref,ralston_struct);

   X_analytical = compute_planetary_motion(t_ref+h_ref, V0, orbit_params);
    
   X_true = [X_true, X_analytical];
   local_error_ralston = [local_error_ralston, norm(XB_ralston - X_analytical)];
   difference = [difference, norm(X_analytical-compute_planetary_motion(t_ref, V0, orbit_params))];

end

figure;
loglog(h_list, local_error_ralston, 'b', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error Ralston Method")
legend("Ralston Method", "Difference")
%% 3/8-Rule - Fourth Order
struct_38.A = [0, 0, 0, 0; 1/3, 0, 0, 0; -1/3, 1, 0, 0; 1, -1, 1, 0];
struct_38.B = [1/8, 3/8, 3/8, 1/8];
struct_38.C = [0; 1/3; 2/3; 1];
h = 0.001;
t0 = 0;          % Start time
tf = 40;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval  

[t_list_38,X_list_38,h_avg_38, num_evals_38] = explicit_RK_fixed_step_integration ...
(my_rate_func,tspan,V0,h,struct_38);

figure()
hold on
plot(X_list_true(:,1), X_list_true(:,2),'b','linewidth',2) 
plot(X_list_38(:,1), X_list_38(:,2), 'r--','linewidth',2)
xlabel("t (sec)")
ylabel("X")
title("3/8-Rule Method Approximated Solution")
legend("3/8-Rule Method", "True Solution")

% Local Error Calculations
h_list = logspace(-5,1,100); 

local_error_eighth = [];
X_true = [];
difference = [];

t_ref = 4.49;
XA = compute_planetary_motion(t_ref, V0, orbit_params);

for i = 1:length(h_list)
   h_ref = h_list(i);
    
   [XB_eighth,~] = explicit_RK_step(my_rate_func,t_ref,XA,h_ref,struct_38);

   X_analytical = compute_planetary_motion(t_ref+h_ref, V0, orbit_params);
    
   X_true = [X_true, X_analytical];
   local_error_eighth = [local_error_eighth, norm(XB_eighth - X_analytical)];
   difference = [difference, norm(X_analytical-compute_planetary_motion(t_ref, V0, orbit_params))];

end


figure;
loglog(h_list, local_error_eighth, 'b', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error 3/8-Rule Method")
legend("3/8-Rule Method", "Difference")

%% LOCAL PLOT
figure;
loglog(h_list, local_error_heun, 'b', 'LineWidth', 2);
hold on
loglog(h_list, local_error_ralston, 'g', 'LineWidth', 2);
hold on
loglog(h_list, local_error_eighth, 'r', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'm', 'LineWidth', 2)
title("Local Error Method Comparison")
xlabel("Step Size")
ylabel("Local Error")
legend("Heun's Method", "Ralston's Method", "3/8-Rule Method", "Difference")

%% explicit_RK fixed_step integrator
%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct)
    % calculate steps and h
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    % define variables
    XA = X0;
    num_evals = 0;
    t_list = linspace(tspan(1),tspan(2),num_steps+1);
 
    X_list = zeros(num_steps+1,length(X0));
    X_list(1,:) = X0';
    %calculate the values until it is just short of the end value
    for i = 1:num_steps
        t = t_list(i);
        [XB, temp_eval] = explicit_RK_step(rate_func_in,t,XA,h_avg,BT_struct);
        num_evals = num_evals + temp_eval;

        X_list(i+1,:)= XB';
        XA = XB;
    end  
end

%% explicit_RK_step
%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    k = zeros(length(XA),length(BT_struct.B));
    for i = 1:length(BT_struct.B)
        k(:,i) = rate_func_in(t+BT_struct.C(i)*h, XA+h*(k*BT_struct.A(i,:)'));
    end
    XB = XA + h*(k*BT_struct.B');
    num_evals = length(BT_struct.B);
end
%% RATE_FUNC01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end
%% ITERATION SOLVER
function [num_steps, h] = iteration_solver(tspan, h_ref)
    t_range = tspan(2)-tspan(1);
    num_steps = t_range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = t_range/num_steps; % Divide range by steps to get real h
end
