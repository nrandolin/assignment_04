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
    r = sqrt(V(1)^2+V(2)^2);
    dVdt(1) = V(3);
    dVdt(2) = V(4);
    dVdt(3) = -orbit_params.m_sun*orbit_params.G*V(1)/(r^3);
    dVdt(4) = -orbit_params.m_sun*orbit_params.G*V(2)/(r^3);
end