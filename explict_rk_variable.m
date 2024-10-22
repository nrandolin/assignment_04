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


% [p_test1,k_test1] = loglog_fit(h_list,local_error_test1);
% [p_test2,k_test2] = loglog_fit(h_list,local_error_test2);
% [p_diff,k_diff] = loglog_fit(h_list,difference);
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