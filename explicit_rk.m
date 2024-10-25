clear;

BT_struct.A = [0, 0; 0.5, 0];
BT_struct.B = [0, 1];
BT_struct.C = [0; 0.5];
t = 2;
XA = solution01(t);
h = 0.1;
t0 = 0;          % Start time
tf = 0.5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval
X0 = 1;      % Initial conditions

[XB, num_evals] = explicit_RK_step(@rate_func01,t,XA,h,BT_struct);

[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,X0,h,BT_struct);

figure;
plot(t_list,X_list);
hold on;
plot(t_list,solution01(t_list))

%% PLANETARY MOTION
% testing compute_planetary motion
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 6;
y0 = 10;
dxdt0 = 0;
dydt0 = 1.5;

V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,40,500);
%V_list = compute_planetary_motion(t_range,V0,orbit_params);
[V_list,E_list,H_list] = compute_planetary_motionVEH(t_range,V0,orbit_params);

figure;
axis equal; axis square;
axis([-20,20,-20,20])
hold on
plot(0,0,'ro','markerfacecolor','r','markersize',5);
plot(V_list(:,1),V_list(:,2),'k');
hold off

figure;
hold on;
axis([0,30,-2,2])
plot(t_range,E_list,"k");
xlabel("time (s)")
ylabel("Mechanical Energy (Joules)")
yyaxis right
ylabel("Angular Momentum (kgm^2/s)")
ylim([-2,2])
plot(t_range,H_list);
title("Conservation of Mechanical Energy and Angular Momentum")
hold off

%% Heun's method
filterparams.min_xval = 0;
filterparams.max_xval = 0.01;

BT_struct.A = [0, 0; 1, 0];
BT_struct.B = [0.5, 0.5];
BT_struct.C = [0; 1];
h = 0.1;
t0 = 0;          % Start time
tf = 5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval  
XA = solution01(t0);


% APPROXIMATION
[t_list_heun,X_list_heun,h_avg_heun,num_evals_heun] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,XA,h,BT_struct);

figure;
plot(t_list_heun,X_list_heun, 'b', 'LineWidth', 2);
hold on;
plot(t_list_heun,solution01(t_list_heun), '--r', 'LineWidth', 2)
title("Heun's Method Approximated Solution")
legend("Heun's Method", "True Solution")

% LOCAL
h_list = logspace(-5,1,100); 

local_error_heun = [];
X_true = [];
difference = [];

t_ref = 4.49;
XA = solution01(t_ref);

for i = 1:length(h_list)
   h_ref = h_list(i);
    
   [XB_heun,~] = explicit_RK_step(@rate_func01,t_ref,XA,h_ref,BT_struct);

   X_analytical = solution01(t_ref+h_ref);
    
   X_true = [X_true, X_analytical];
   local_error_heun = [local_error_heun, norm(XB_heun - X_analytical)];
   difference = [difference, norm(X_analytical-solution01(t_ref))];

end

figure;
[p_heun,k_heun] = loglog_fit(h_list,local_error_heun)
[p_diff,k_diff] = loglog_fit(h_list,difference)
loglog(h_list, local_error_heun, 'b', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error Heun's Method")
legend("Heun Method's", "Difference")

% GLOBAL ERROR
% h_list = logspace(-3,1,50);  % Time step sizes
% global_error_heun = [];
% rate_function_calls_heun = [];
% X_true = [];
% results_matrix = zeros(length(h_list), 6);
% 
% for i = 1:length(h_list)
%    h_ref = h_list(i);
%    [t_list_heun,X_list_heun,h_avg_heun,num_evals_heun] = explicit_RK_fixed_step_integration ...
%    (@rate_func01,tspan,XA,h,BT_struct);
%   
%    X_numerical = X_list_heun(end,:)'; 
%    X_analytical = solution01(tf);
% 
% 
%    global_error_heun = [global_error_heun, norm(X_numerical - X_analytical)];
%    rate_function_calls_heun = [rate_function_calls_heun, num_evals_heun];
% 
% end

%[p_heun,k_heun] = loglog_fit(rate_function_calls_heun,global_error_heun, filterparams);
%[p_heun_step,k_heun_step] = loglog_fit(h_list,global_error_heun, filterparams);

% figure(3);
% loglog(h_list, num_evals_heun);
% loglog(rate_function_calls_heun, k_heun*rate_function_calls_heun.^p_heun, 'r--', 'LineWidth', 1.5);

%% Ralston's Third-Order Method
filterparams.min_xval = 0;
filterparams.max_xval = 0.01;
BT_struct.A = [0, 0, 0; 0.5, 0, 0; 0, 0.75, 0];
BT_struct.B = [2/9, 1/3, 4/9];
BT_struct.C = [0; 0/5; 0.75];
t0 = 0;          % Start time
tf = 5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval  
XA = solution01(t0);
h = 0.1;

[t_list_ralston,X_list_ralston,h_avg_ralston,num_evals_ralston] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,XA,h,BT_struct);

figure;
plot(t_list_ralston,X_list_ralston, 'b', 'LineWidth', 2);
hold on;
plot(t_list_ralston,solution01(t_list_ralston), '--r', 'LineWidth', 2)
title("Ralston's Third-Order Method Approximated Solution")
legend("Ralston's Third-Order Method", "True Solution")

% LOCAL
h_list = logspace(-5,1,100); 

local_error_ralston = [];
X_true = [];
difference = [];

t_ref = 4.49;
XA = solution01(t_ref);

for i = 1:length(h_list)
   h_ref = h_list(i);
    
   [XB_ralston,~] = explicit_RK_step(@rate_func01,t_ref,XA,h_ref,BT_struct);

   X_analytical = solution01(t_ref+h_ref);
    
   X_true = [X_true, X_analytical];
   local_error_ralston = [local_error_ralston, norm(XB_ralston - X_analytical)];
   difference = [difference, norm(X_analytical-solution01(t_ref))];

end

figure;
[p_ralston,k_ralston] = loglog_fit(h_list,local_error_ralston);
[p_diff,k_diff] = loglog_fit(h_list,difference);
loglog(h_list, local_error_ralston, 'b', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error Ralston Method")
legend("Ralston Method", "Difference")

% GLOBAL
% 
% h_list = logspace(-3,1,50);  % Time step sizes
% global_error_ralston = [];
% rate_function_calls_ralston = [];
% X_true = [];
% results_matrix = zeros(length(h_list), 6);
% 
% for i = 1:length(h_list)
%    h_ref = h_list(i);
%    [t_list_ralston,X_list_ralston,h_avg_ralston,num_evals_ralston] = explicit_RK_fixed_step_integration ...
%    (@rate_func01,tspan,XA,h,BT_struct);
%   
%    X_numerical = X_list_ralston(end,:)'; 
%    X_analytical = solution01(tf);
% 
% 
%    global_error_ralston = [global_error_ralston, norm(X_numerical - X_analytical)];
%    rate_function_calls_ralston = [rate_function_calls_ralston, num_evals_heun];
% 
% end
% 
% [p_g_ralston,k_g_ralston] = loglog_fit(rate_function_calls_ralston,global_error_ralston, filterparams);
% [p_ralston_step,k_ralston_step] = loglog_fit(h_list,global_error_ralston, filterparams);
% 
% figure(3);
% loglog(h_list, num_evals_ralston);
% loglog(rate_function_calls_ralston, k_ralston_step*rate_function_calls_ralston.^p_ralston_step, 'r--', 'LineWidth', 1.5);

%% 3/8-Rule Fourth-Order Method
BT_struct.A = [0, 0, 0, 0; 1/3, 0, 0, 0; -1/3, 1, 0, 0; 1, -1, 1, 0];
BT_struct.B = [1/8, 3/8, 3/8, 1/8];
BT_struct.C = [0; 1/3; 2/3; 1];
t0 = 0;          % Start time
tf = 5;        % End time (replace with the desired final time)
tspan = [t0, tf];  % Full integration interval  
XA = solution01(t0);
h = 0.1;

[t_list_eighth,X_list_eighth,h_avg_eighth,num_evals_eighth] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,XA,h,BT_struct);

figure;
plot(t_list_eighth,X_list_eighth, 'b', 'LineWidth', 2);
hold on;
plot(t_list_eighth,solution01(t_list_eighth), '--r', 'LineWidth', 2)
title("3/8-Rule Fourth-Order Approximated Solution")
legend("3/8-Rule Fourth-Order", "True Solution")

% LOCAL
h_list = logspace(-5,1,100); 

local_error_eighth = [];
X_true = [];
difference = [];

t_ref = 4.49;
XA = solution01(t_ref);

for i = 1:length(h_list)
   h_ref = h_list(i);
    
   [XB_eighth,~] = explicit_RK_step(@rate_func01,t_ref,XA,h_ref,BT_struct);

   X_analytical = solution01(t_ref+h_ref);
    
   X_true = [X_true, X_analytical];
   local_error_eighth = [local_error_eighth, norm(XB_eighth - X_analytical)];
   difference = [difference, norm(X_analytical-solution01(t_ref))];

end


[p_eighth,k_eighth] = loglog_fit(h_list,local_error_eighth);
[p_diff,k_diff] = loglog_fit(h_list,difference);
figure;
loglog(h_list, local_error_eighth, 'b', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error 3/8-Rule Method")
legend("3/8-Rule Method", "Difference")

%% FULL LOCAL PLOT
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

        X_list(i+1,:)= XB;
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
    range = tspan(2)-tspan(1);
    num_steps = range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = range/num_steps; % Divide range by steps to get real h
end
