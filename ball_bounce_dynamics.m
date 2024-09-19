function dXdt = ball_bounce_dynamics(t,X,params)
% second-order ODEs reduced to first order
% X == state vector
%    y = X(1) == ball vertical position
%    y_dot = X(2) == ball vertical velocity
% params is a struct array containing static model parameters like mass and
% contact parameters
y = X(1);
y_dot = X(2);

% For this bouncing model, we need to calculate a contact force using hunt-crossley model (only engages
% when in contact with ground, delta > 0)
k = params.k;
a = params.a;
radius = params.radius;
delta = -(y-radius); 
delta_dot = -y_dot;
if delta > 0
    F = k*delta*(3/2)*(1+a*delta_dot);
else
    F = 0;
end

if F < 0 % enforce unilateral contact in the case of positive delta but large negative delta_dot
    F = 0;
end

% Here, we finally evaluate our differential equation at the input t and X
m = params.m;
g = params.g;
dXdt(1,1) = y_dot; % time derivative of y = y_dot;
dXdt(2,1) = -g + F/m; % time derivative of y_dot is N2L (our 2nd order ODE)
% NOTE: always return dXdt as a column vector!!

end