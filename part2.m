
%  Numerical integration of equations of motion for
%  golf ball trajectory (derived in class).
%
%  x(1)=Vx, x(2)=Vy, x(3)=Vz, x(4)=X, x(5)=Y, x(6)=Z, x(7)=omega
%
%  Updated:  30-Sept-2019  by  J. McPhee
%

global radius mass rho area inertia grav tx ty tz

pie = 3.1415926535;
grav = 32.17;
radius = (1.68/2)/12;    % diameter of 1.68 inches
mass = (1.62/16)/grav;   % mass of 1.62 ounces
rho_arr = 0.0023769*[0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5];         % density of air
area = pie*radius*radius;
inertia = 0.4*mass*radius*radius;   % inertia of a sphere


params.static_Cd = Cd; 
params.static_Cl = Cl;
params.static_Cm = Cm;
params.radius = radius;
params.mass = mass;
params.rho = rho;
params.area = area;
params.inertia = inertia;
params.grav = grav;

%  Launch conditions
%Shot 10 is used for this example
carry_values = zeros([1 length(rho_arr)]);
for i=1:length(rho_arr)
    rho = rho_arr(i);
    v0 = 163; % in mph
    elev = 14.4; % in deg
    azim = -0.9;  % in deg, about Y (i.e. going left)
    back = 4184; % in rpm, about global Z
    side = 51;    % in rpm, about global Y
    rifle = 0;   % in rpm, about global X

    vx = (v0*cos(elev*pie/180)*cos(azim*pie/180))*88/60;  % convert mph to ft/s
    vy = v0*sin(elev*pie/180)*88/60;
    vz = -(v0*cos(elev*pie/180)*sin(azim*pie/180))*88/60;
    wx = rifle*pie/30;   % convert rpm to rad/s
    wy = side*pie/30;
    wz = back*pie/30;

    omega = sqrt(wx*wx + wy*wy + wz*wz);
    tx = wx/omega;
    ty = wy/omega;
    tz = wz/omega;
    params.tx = tx; % add to param struct
    params.ty = ty;
    params.tz = tz;
    
    t0 = 0;
    tf = 10;
    x0 = [vx, vy, vz, 0, 0, 0, omega]';   % launch conditions
    %x0 = [0, 6, 0, vx, vy, vz, omega]';

    dynamic_eqn = @(t,x) golf_eqns(x, params);
    options = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t,x] = ode45(dynamic_eqn, [t0,tf], x0, options);

    if x(end,5) > 0
        % ball still in the air
        disp('Ball still in the air, consider changing tf')
        X = x(1:end,4)/3;
        Y = x(1:end,5)/3;
        Z = x(1:end,6)/3;
    else
        ground = find(x(:,5) < 0, 1);
        x_ground = interp1(x(ground-1:ground, 5), x(ground-1:ground,:), 0);
        x_ground(5) = 0;
        x = x(1:ground, :);
        x(end,:) = x_ground;
        X = x(1:end,4)/3;
        Y = x(1:end,5)/3;
        Z = x(1:end,6)/3;
    end

    carry = X(end); % X driving distance in yards
    carry_values(i) = carry;
   
end

figure
plot(rho_arr, carry_values)
xlabel('Air Density (slug/ft^3)')
ylabel('Carry Distance (yards)')


