%fits the relationship y = k*x^p to (x,y) data
%INPUTS:
%x_regression: the x data (a row or column vector)
%y_regression: the y data (a row or column vector)
%filter_params (optional): a struct with parameters to set the range of
%data points that should be included in the regression
%you can set as many or as few of these parameters as you'd like
%   filter_params.min_index
%   filter_params.max_index
%   filter_params.min_xval
%   filter_params.max_xval
%   filter_params.min_yval
%   filter_params.max_yval
%OUTPUTS:
%[p,k]: the regressed values for relationship y = k*x^p
function [p,k] = loglog_fit(x_regression,y_regression,varargin)
    
    %convert x_regression to a column vector if it's a row vector
    if size(x_regression,1)==1
        x_regression = abs(x_regression)';
    end

    %convert y_regression to a column vector if it's a row vector
    if size(y_regression,1)==1
        y_regression = abs(y_regression)';
    end

    %if filter_params has been provided, then filter the data points
    if nargin==3
        filter_params = varargin{1};
        
        num_points = length(x_regression);
        indices = 1:num_points;

        filter_bool = ones(num_points,1);

        if isfield(filter_params,'min_index')
            filter_bool = filter_bool & indices>=filter_params.min_index;
        end

        if isfield(filter_params,'max_index')
            filter_bool = filter_bool & indices<=filter_params.max_index;
        end

        if isfield(filter_params,'min_xval')
            filter_bool = filter_bool & x_regression>=filter_params.min_xval;
        end

        if isfield(filter_params,'max_xval')
            filter_bool = filter_bool & x_regression<=filter_params.max_xval;
        end

        if isfield(filter_params,'min_yval')
            filter_bool = filter_bool & y_regression>=filter_params.min_yval;
        end

        if isfield(filter_params,'max_yval')
            filter_bool = filter_bool & y_regression<=filter_params.max_yval;
        end

        x_regression = x_regression(filter_bool);
        y_regression = y_regression(filter_bool);
    end

    %compute the logs of x_regression and y_regression
    Y = log(y_regression);
    X1 = log(x_regression);

    %set up the regression
    X2 = ones(length(X1),1);

    %run the regression
    coeff_vec = regress(Y,[X1,X2]);
    
    %pull out the coefficients from the fit
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end