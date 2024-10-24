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
h_ref = 0.01;
[XB, num_evals, h_next, redo] = explicit_RK_variable_step...
(@rate_func01,t,XA,h,DormandPrince,p,error_desired);
[t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
(@rate_func01,tspan,X0,h_ref,DormandPrince,p,error_desired)
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
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)
    num_evals = 0;
    t = tspan(1);
    h_test = [h_ref, h_ref];
    h_list = [];
    X_list = X0;
    t_list = t; % change num_steps
    XA = X0;
    redo = 1;
    % probably need a for loop for iterating through time, not sure how to
    % do this with varying step size
    while t ~= tspan(2)
        while redo == 1
            [XB, num_evals_temp, h_next, redo] = explicit_RK_variable_step...
            (rate_func_in,t,XA,h_test(end),BT_struct,p,error_desired);
            h_test = [h_test, h_next]; %add the next predicted h to the loop
        end
        h = h_test(end-1);  % the actual h used was one less than the one now
        t_temp = t+h;
        % end early?
        if t_temp > tspan(2)
            h_final = tspan(2)-t_list(end);
            h_test = [h_final, h_final];
            continue
        end
        t = t_temp;
        h_list = [h_list, h]; % create a list of used h values
        t_list = [t_list, t]; %add timestep to list
        X_list = [X_list, XB];% update evaluation
        redo = 1;
        num_evals = num_evals + num_evals_temp;
    end

    h_avg = mean(h_list);

%     % calculate steps and h
%     [num_steps, h_avg] = iteration_solver(tspan, h_ref);
%     % define variables
%     XA = X0;
%     num_evals = 0;
%     t_list = linspace(tspan(1),tspan(2),num_steps+1);
%  
%     X_list = zeros(num_steps+1,length(X0));
%     X_list(1,:) = X0';
%     %calculate the values until it is just short of the end value
%     for i = 1:num_steps
%         t = t_list(i);
%         [XB, temp_eval] = explicit_RK_step(rate_func_in,t,XA,h_avg,BT_struct);
%         num_evals = num_evals + temp_eval;
% 
%         X_list(i+1,:)= XB;
%         XA = XB;
%     end  
end

%% RK Variable Step
function [XB, num_evals, h_next, redo] = explicit_RK_variable_step...
(rate_func_in,t,XA,h,BT_struct,p,error_desired)
    alpha = 4; % btwn 1.5 and 10, inclusive
    [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct); %run 1 step of the solver (on original ts)
    h_next = h*min(0.9*(error_desired/abs(XB1-XB2))^(1/p),alpha); % calculate h_next
    XB = XB1;
    estimated_error = abs(XB1 - XB2);% calculate error
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