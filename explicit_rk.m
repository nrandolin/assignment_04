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

figure(1);
plot(t_list,X_list);
hold on;
plot(t_list,solution01(t_list))

% testing compute_planetary motion
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 8;
y0 = 0;
dxdt0 = 0;
dydt0 = 1.5;

V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,30,100);
V_list = compute_planetary_motion(t_range,V0,orbit_params);


axis equal; axis square;
axis([-20,20,-20,20])
hold on
plot(0,0,'ro','markerfacecolor','r','markersize',5);
plot(V_list(:,1),V_list(:,2),'k');

%% Heun's method
BT_struct.A = [0, 0; 1, 0];
BT_struct.B = [0.5, 0.5];
BT_struct.C = [0; 1];
t = 2;
XA = solution01(t);
h = 0.1;

[t_list_heun,X_list_heun,h_avg_heun,num_evals_heun] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,X0,h,BT_struct);

figure(1);
plot(t_list_heun,X_list_heun, 'b');
hold on;
plot(t_list_heun,solution01(t_list_heun), 'r')

%% Ralston's Third-Order Method
BT_struct.A = [0, 0, 0; 0.5, 0, 0; 0, 0.75, 0];
BT_struct.B = [2/9, 1/3, 4/9];
BT_struct.C = [0; 0/5; 0.75];
t = 2;
XA = solution01(t);
h = 0.1;

[t_list_ralston,X_list_ralston,h_avg_ralston,num_evals_ralston] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,X0,h,BT_struct);

figure(1);
plot(t_list_ralston,X_list_ralston, 'b');
hold on;
plot(t_list_ralston,solution01(t_list_ralston), 'r')

%% 3/8-Rule Fourth-Order Method
BT_struct.A = [0, 0, 0, 0; 1/3, 0, 0, 0; -1/3, 1, 0, 0; 1, -1, 1, 0];
BT_struct.B = [1/8, 3/8, 3/8, 1/8];
BT_struct.C = [0; 1/3; 2/3; 1]
t = 2;
XA = solution01(t);
h = 0.1;

[t_list_eighth,X_list_eighth,h_avg_eighth,num_evals_eighth] = explicit_RK_fixed_step_integration ...
(@rate_func01,tspan,X0,h,BT_struct);

figure(1);
plot(t_list_eighth,X_list_eighth, 'b');
hold on;
plot(t_list_eighth,solution01(t_list_eighth), 'r')

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
    k = zeros(length(BT_struct.B), 1);
    num_loop = 0;
    for i = 1:length(BT_struct.B)
        for j = 1:length(BT_struct.B)-1
            k(i) = k(i) + rate_func_in(t, t + BT_struct.C(i)*h*XA + h*BT_struct.A(i,j)*k(j));
            num_loop = num_loop+1;
        end 
    end
    func = @(XB) XA + h*BT_struct.B*k;
    [XB, num_evals] = multi_newton_solver(func, XA, true);
    num_evals = num_evals + num_loop;
end

%% MULTI NEWTON SOLVER
function [x_next, num_evals] = multi_newton_solver(fun,x_guess,varargin)
%true if supposed to use analytical jacobian, false otherwise
use_analytical_jacobian = nargin==3 && varargin{1}(1);

A_thresh = 10e-14;
B_thresh = 10e-14;
num_evals = 0;

    if use_analytical_jacobian == true
    f_val = fun(x_guess);
    num_evals = 1;
    J = approximate_jacobian(fun,x_guess);
    x_next = x_guess;
        while max(abs(J\f_val))> A_thresh && max(abs(f_val))>B_thresh && abs(det(J*J')) > 1e-14
            %calculate next x
            x_next = x_next - J\f_val;
            %establish next y value
            f_val = fun(x_next);
            J = approximate_jacobian(fun,x_next);  
            num_evals = num_evals+1;
        end
    end

%Loop through until x is small
    if use_analytical_jacobian == false
    [f_val, J] = fun(x_guess);
    x_next = x_guess;
        while max(abs(J\f_val))> A_thresh && max(abs(f_val))>B_thresh && abs(det(J*J')) > 1e-14
            %calculate next x
            x_next = x_next - J\f_val;
            %establish next y value
            [f_val, J]  = fun(x_next);        
        end
    end

end
%% APPROX JACOBIAN
function J = approximate_jacobian(fun,x)
    % Set initial variables
    ej = zeros(length(x), 1); %variable to store vector of multiplyers
    h = 1e-6;
    J = zeros(length(fun(x)), length(x));
    for i = 1:size(J, 2)
        ej(i) = 1;
        % calculate the partial derivative 
        J(:, i) = (fun(x+h*ej) - fun(x-h*ej))/(2*h);
        ej(i) = 0;
    end 
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
%% Universal Gravitation rate func
%Rate function describing Newton’s law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
    dVdt = zeros(4,1);
    dVdt(1) = V(3, t);
    dVdt(2) = V(4, t);
    dVdt(3) = -orbit_params.m_sun*orbit_params.G*V(1, t)/(abs(V(1, t))^3);
    dVdt(4) = -orbit_params.m_sun*orbit_params.G*V(2, t)/(abs(V(2, t))^3);
end
