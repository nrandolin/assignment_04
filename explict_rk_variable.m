DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
DormandPrince.A = [0,0,0,0,0,0,0;
1/5, 0, 0, 0,0,0,0;...
3/40, 9/40, 0, 0, 0, 0,0;...
44/45, -56/15, 32/9, 0, 0, 0,0;...
19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

t = 2;
XA = solution01(t);
h = 0.03;
[XB1, XB2, num_evals] = RK_step_embedded(@rate_func01,t,XA,h,DormandPrince);
%% testing variable step
p = 3;
error_desired = 0.05;
tspan = [0, 2];
X0 = XA;
h_ref = 0.05;
[XB, num_evals, h_next, redo] = explicit_RK_variable_step...
(@rate_func01,t,XA,h,DormandPrince,p,error_desired);
[t_list,X_list,h_avg, num_evals, percent_failed] = explicit_RK_variable_step_integration ...
(@rate_func01,tspan,X0,h_ref,DormandPrince,p,error_desired);

%% testing integrator
desired_error_range = linspace(0.0001,0.001,11);
x0 = 6;
y0 = 10;
dxdt0 = 0;
dydt0 = 1.5;
h=0.03;
p=3;
h_ref = 0.03;
V0 = [x0;y0;dxdt0;dydt0];
tspan = [0, 30];
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
my_rate_func = @(t_in,V_in) gravity_rate_func(t_in,V_in,orbit_params);
h_list = linspace(0.001, 0.15, 100);
error_list = logspace(-12, -1, 100);
h_avg_list = [];
global_error_list = [];
percent_failed_list = [];
num_evals_list = [];
for i = 1:length(error_list)
    error = error_list(i);
    [t_list,X_list,h_avg, num_evals, percent_failed] = explicit_RK_variable_step_integration ...
    (my_rate_func,tspan,V0,h_ref,DormandPrince,p,error);
    X_numerical = X_list(:, end);
    X_analytical = compute_planetary_motion(tspan(2),V0,orbit_params);
    global_error = norm(X_numerical - X_analytical);
    global_error_list = [global_error_list, global_error];
    h_avg_list = [h_avg_list, h_avg];
    percent_failed_list = [percent_failed_list, percent_failed];
    num_evals_list = [num_evals_list, num_evals];
end
V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,40,500);
X_list_true = compute_planetary_motion(t_range,V0,orbit_params);

h_list = linspace(0.001, 0.12, 100);
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
%[XB, num_evals] = explicit_RK_step(@rate_func01,t,XA,h,BT_struct);
global_error_list_fixed = [];
h_avg_list_fixed = [];
num_evals_list_fixed = [];
for i = 1:length(h_list)
    h = h_list(i);
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
    (my_rate_func,tspan,V0,h,DormandPrince);
    X_numerical = X_list(end, :)';
    X_analytical = compute_planetary_motion(tspan(2),V0,orbit_params);
    global_error = norm(X_numerical - X_analytical);
    global_error_list_fixed = [global_error_list_fixed, global_error];
    h_avg_list_fixed = [h_avg_list_fixed, h_avg];
    num_evals_list_fixed = [num_evals_list_fixed, num_evals];
end


% global error vs. avg step size
figure()
loglog(h_avg_list_fixed, global_error_list_fixed, '.b', 'MarkerSize', 8)
hold on
loglog(h_avg_list, global_error_list, '.r', 'MarkerSize', 8)
title("Global Error vs. Average Step Size")
xlabel("Average Step Size")
ylabel("Global Error")
legend("Fixed Step", "Variable Step")

% global error vs. num evals
figure()
loglog(num_evals_list_fixed, global_error_list_fixed, '.b', 'MarkerSize', 8)
hold on
loglog(num_evals_list, global_error_list, '.r', 'MarkerSize', 8)
title("Global Error vs. Number of Evaluations")
xlabel("Number of Evaluations")
ylabel("Global Error")
legend("Fixed Step", "Variable Step")

figure()
plot(h_avg_list, percent_failed_list, '.m', 'MarkerSize', 8)
set(gca, 'XScale', 'log')
title("Percent Failed vs. Average Step Size")
xlabel("Average Step Size")
ylabel("Percent Failed")


% plot planet trajectory
figure()
hold on
plot(X_list_true(:,1), X_list_true(:,2),'b','linewidth',2) 
plot(X_list(:,1), X_list(:,2), 'r--','linewidth',2)
title("Approximated vs. True Planet Trajectory")
legend("True Solution", "Approximated Solution")

%% Velocity and Position
tspan = [0, 40];
DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];

[t_list,X_list,h_avg, num_evals, percent_failed] = explicit_RK_variable_step_integration ...
        (my_rate_func,tspan,V0,h_ref,DormandPrince,p,.0001);

% position vs. time
figure()
plot(t_list, X_list(1,:), 'bo-','markerfacecolor','k','markeredgecolor','k','markersize',2)
hold on
plot(t_list, X_list(2,:), 'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2)
title('Position vs. Time')
legend('x', 'y')
xlabel('Time (s)')
ylabel('Position (m)')

% velocity vs. time
figure()
plot(t_list, X_list(3,:), 'bo-','markerfacecolor','k','markeredgecolor','k','markersize',2)
hold on
plot(t_list, X_list(4,:), 'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2)
title('Velocity vs. Time')
legend('x', 'y')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

%% step size vs distance
tspan = [0, 40];
h_ref = 10;

[t_list,X_list,h_avg, num_evals, percent_failed] = explicit_RK_variable_step_integration ...
        (my_rate_func,tspan,V0,h_ref,DormandPrince,p,.0001);
X_list
h_list = [diff(t_list)];
radius = sqrt(X_list(1,:).^2 + X_list(2,:).^2);
radius(1) = [];

% figure()
% 
% for i = 1:length(radius)
%     xlim([0 20])
%     ylim([0 5])
%     plot(radius(i), h_list(i), 'ob');
%     hold on;
%     pause(0.25);
% end

figure()
loglog(radius, h_list, 'o')
title("Step Size vs. Distance from Planet to Sun")
ylabel("Step Size (s)")
xlabel("Distance (m)")



%%
filterparams.min_xval = 0;
filterparams.max_xval = 0.1;
h_list = logspace(-3,1,50);  % Time step sizes
global_error = [];
rate_function_calls = [];
X_true = [];
tspan = [0, 40];
results_matrix = zeros(length(h_list), 6);
tf = tspan(end);
X_list = [];
for i = 1:length(h_list)
   h_ref = h_list(i);
   [t_list,X_list_temp,h_avg,num_evals_temp] = explicit_RK_variable_step_integration ...
   (my_rate_func,tspan,V0,h_ref,DormandPrince,p,.0001);
   X_numerical = X_list_temp(:, end)';
   X_analytical = compute_planetary_motion(tf,V0,orbit_params);
   global_error = [global_error, norm(X_numerical - X_analytical)];
   rate_function_calls = [rate_function_calls, num_evals_temp];


end
%end



%p_heun_step,k_heun_step] = loglog_fit(h_list,global_error, filterparams);
%[p_heun,k_heun] = loglog_fit(rate_function_calls,global_error, filterparams);

% figure(3);
% figure();
% loglog(h_list, global_error, 'b');
%loglog(rate_function_calls, k_heun*rate_function_calls.^p_heun, 'r--', 'LineWidth', 1.5);

%% LOCAL ERROR
h_list = logspace(-5,1,100); 

local_error_test1 = [];
local_error_test2 = [];
X_true = [];
difference = [];
minus = [];

t_ref = 4.49;
XA = solution01(t_ref);

for i = 1:length(h_list)
   h_ref = h_list(i);
    
   [XB1, XB2, ~] = RK_step_embedded(@rate_func01,t_ref,XA,h_ref,DormandPrince);

   X_analytical = solution01(t_ref+h_ref);
    
   X_true = [X_true, X_analytical];
   local_error_test1 = [local_error_test1, norm(XB1 - X_analytical)];
   local_error_test2 = [local_error_test2, norm(XB2 - X_analytical)];
   difference = [difference, norm(X_analytical-solution01(t_ref))];
   minus = [minus, abs(XB1-XB2)];

end

figure;
loglog(h_list, local_error_test1, 'b', 'LineWidth', 2);
hold on
loglog(h_list, local_error_test2, 'g', 'LineWidth', 2);
hold on
loglog(h_list, minus, 'm', 'LineWidth', 2);
hold on
loglog(h_list, difference, 'r', 'LineWidth', 2)
title("Local Error Test")
legend("XB1", "XB2", "|XB1-XB2|","Difference")
ylabel("Local Error")
xlabel("Step Size")

figure;
loglog(minus, local_error_test1, 'b', 'LineWidth', 2);
hold on
loglog(minus, local_error_test2, 'g', 'LineWidth', 2);
title("Local Error vs. |XB1-XB2|")
legend("XB1", "XB2")
ylabel("Local Error")
xlabel("|XB1-XB2|")
%% RK Variable Step Integration
%Runs numerical integration arbitrary RK method using variable time steps
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
%p: how error scales with step size (error = k*h^p)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, percent_failed] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)
    num_evals = 0;
    t = tspan(1);
    tf = tspan(2);
%     h_test = [h_ref, h_ref];
%     h_list = [];
    X_list = X0;
    t_list = t; % change num_steps
    XA = X0;
%     redo = 1;

    h = h_ref;

    num_failed_steps = 0;

    num_attempted_steps = 0;
    while t<tf
        num_attempted_steps = num_attempted_steps+1;
        t_next = t+h;

        if t_next>tf
            h= tf-t;
            t_next = tf;
        end

        [XB, num_evals_temp, h_next, redo] = explicit_RK_variable_step...
                (rate_func_in,t,XA,h,BT_struct,p,error_desired);

        num_evals = num_evals+num_evals_temp;
        h = h_next;
        
        if ~redo
            XA = XB;
            t = t_next;
            X_list(:,end+1) = XA;
            t_list(end+1) = t;
        else
            num_failed_steps = num_failed_steps+1;
        end

    end

    h_avg = (tspan(2)-tspan(1))/(length(t_list)-1);
    percent_failed = num_failed_steps/num_attempted_steps;


    % probably need a for loop for iterating through time, not sure how to
    % do this with varying step size
%     while t ~= tspan(2)
%         %for i = 1:length(XA)
%             while redo == 1
%                 num_failed_steps = num_failed_steps + 1;
%                 [XB, num_evals_temp, h_next, redo] = explicit_RK_variable_step...
%                 (rate_func_in,t,XA,h_test(end),BT_struct,p,error_desired);
%                 h_test = [h_test, h_next]; %add the next predicted h to the loop
%             end
%         %end
%         % Not sure if we want this because it seems to make failure = 0%
% %         if num_failed_steps > 0
% %             num_failed_steps = num_failed_steps - 1; 
% %         end
%         h = h_test(end-1);  % the actual h used was one less than the one now
%         t_temp = t+h;
%         % end early?
%         if t_temp > tspan(2)
%             h_final = tspan(2)-t_list(end);
%             h_test = [h_final, h_final];
%             continue
%         end
%         t = t_temp;
%         h_list = [h_list, h]; % create a list of used h values
%         t_list = [t_list, t]; %add timestep to list
%         X_list = [X_list, XB];% update evaluation
%         
%         redo = 1;
%         num_evals = num_evals + num_evals_temp;
%     end
%     percent_failed = num_failed_steps/num_evals;
%     h_avg = mean(h_list);

end

%% RK Variable Step
function [XB, num_evals, h_next, redo] = explicit_RK_variable_step...
(rate_func_in,t,XA,h,BT_struct,p,error_desired)
    alpha = 4; % btwn 1.5 and 10, inclusive
    [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct); %run 1 step of the solver (on original ts)
    h_next = h*min(0.9*(error_desired/norm(XB1-XB2))^(1/p),alpha); % calculate h_next
    XB = XB1;
    estimated_error = norm(XB1 - XB2); % calculate error
    redo = error_desired<estimated_error;
end
%% RK_step_embedded
%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
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
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    k = zeros(length(XA),length(BT_struct.B));
    for i = 1:length(BT_struct.B)
        k(:,i) = rate_func_in(t+BT_struct.C(i)*h, XA+h*(k*BT_struct.A(i,:)'));
    end
    XB1 = XA + h*(k*BT_struct.B(1,:)');
    XB2 = XA + h*(k*BT_struct.B(2,:)');
    num_evals = length(BT_struct.B);
end

%% RATE_FUNC01
function dXdt = rate_func01(t,X)
dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
X = cos(t);
end
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
%% ITERATION SOLVER
function [num_steps, h] = iteration_solver(tspan, h_ref)
    t_range = tspan(2)-tspan(1);
    num_steps = t_range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = t_range/num_steps; % Divide range by steps to get real h
end