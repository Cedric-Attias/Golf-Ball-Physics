function xdot = golf_eqns_2(t,x)

%  Numerical integration of equations of
%  motion for golf ball using Quintavalla model
%
%  x(1)=Vx, x(2)=Vy, x(3)=Vz, x(4)=X, x(5)=Y, x(6)=Z, x(7)=omega

xdot = zeros(12,1);   % pre-define to make Matlab run faster

global radius mass rho area inertia grav tx ty tz

speed = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3));  % total speed
omega = x(7);                                 % total spin
spinratio = radius*omega/speed;
Q = rho*speed*speed*area/2;

% Cd = 0.28;   % 0.2 is high-speed value from Smith 2018
% Cl = 0.25;     % 0.2 is average value from Smith 2018
% Cm = 0.1;

% Cd = 0.171 + 0.62*spinratio;   % from Steve Quintavalla's paper.
% Cl = 0.083 + 0.885*spinratio;
% Cm = 0.0125*spinratio;

Cd = x(8) + x(9)*spinratio;   % from Steve Quintavalla's paper.
Cl = x(10) + x(11)*spinratio;
Cm = x(12)*spinratio;

xdot(1) = (-Q*Cd*x(1)/speed + Q*Cl*(ty*x(3)-tz*x(2))/speed)/mass ;

xdot(2) = (-Q*Cd*x(2)/speed + Q*Cl*(tz*x(1)-tx*x(3))/speed)/mass - grav;

xdot(3) = (-Q*Cd*x(3)/speed + Q*Cl*(tx*x(2)-ty*x(1))/speed)/mass;

xdot(4) = x(1);
xdot(5) = x(2);
xdot(6) = x(3);

xdot(7) = -Q*Cm*radius*2/inertia;

end


