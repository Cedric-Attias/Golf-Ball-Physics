speed_real = [162 170 171 171 167 161 116 122 167 163];
mean_speed = sum(speed_real)/10;

param_start = [0 1000];
lb = [0, 1000];
ub = [45, 4000];
param_opt = fmincon(@(p) opt(mean_speed, p(1), p(2)) ,param_start,[],[],[],[],lb,ub);

function opt_func = opt(speed, elevation, backspin)

    global radius mass rho area inertia grav tx ty tz

    pie = 3.1415926535;
    grav = 32.17;
    radius = (1.68/2)/12;    % diameter of 1.68 inches
    mass = (1.62/16)/grav;   % mass of 1.62 ounces
    rho = 0.0023769;         % density of air
    area = pie*radius*radius;
    inertia = 0.4*mass*radius*radius;   % inertia of a sphere

    %  Launch conditions
    v0 = speed; % in mph
    elev = elevation; % in deg
    azim = 0;  % in deg, about Y (i.e. going left)
    back = backspin; % in rpm, about global Z
    side = 0;    % in rpm, about global Y
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

    t0 = 0;
    tf = 10;
    x0 = [vx, vy, vz, 0, 0, 0, omega]';   % launch conditions

    options = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t,x] = ode45('golf_eqns_3', [t0,tf], x0, options);

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

    carry = X(end) % X driving distance in yards
    opt_func = -carry;
    
end