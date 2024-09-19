function xdot = golf_eqns_3(t,x)

%  Numerical integration of equations of
%  motion for golf ball using Quintavalla model
%
%  x(1)=Vx, x(2)=Vy, x(3)=Vz, x(4)=X, x(5)=Y, x(6)=Z, x(7)=omega

xdot = zeros(7,1);   % pre-define to make Matlab run faster

global radius mass rho area inertia grav tx ty tz

speed = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3));  % total speed
omega = x(7);                                 % total spin
spinratio = radius*omega/speed;
Q = rho*speed*speed*area/2;

a = 0.171;
b = 0.62;
c = 0.083;
d = 0.885;
e = 0.0125;

Cd = a + b*spinratio;   % from Steve Quintavalla's paper.
Cl = c + d*spinratio;
Cm = e*spinratio;

xdot(1) = (-Q*Cd*x(1)/speed + Q*Cl*(ty*x(3)-tz*x(2))/speed)/mass ;

xdot(2) = (-Q*Cd*x(2)/speed + Q*Cl*(tz*x(1)-tx*x(3))/speed)/mass - grav;

xdot(3) = (-Q*Cd*x(3)/speed + Q*Cl*(tx*x(2)-ty*x(1))/speed)/mass;

xdot(4) = x(1);
xdot(5) = x(2);
xdot(6) = x(3);

xdot(7) = -Q*Cm*radius*2/inertia;

end

