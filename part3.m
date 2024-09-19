
carry_real = [252 286 277 274 283 268 184 183 279 264];
offline_real = [45 -32 -11 -12 -1 -13 -17 -15 -53 -1];

param_start = [0 0 0 0 0];
lb = [0, 0, 0, 0, 0];
ub = [1, 1, 1, 1, 1];
param_opt = fmincon(@(p) LPE(p(1), p(2), p(3), p(4), p(5), carry_real, offline_real) ,param_start,[],[],[],[],lb,ub);

function LPE_func = LPE(a,b,c,d,e,carry_r, offline_r)

global radius mass rho area inertia grav tx ty tz

% a = param_arr(1);
% b = param_arr(2);
% c = param_arr(3);
% d = param_arr(4);
% e = param_arr(5);

pie = 3.1415926535;
grav = 32.17;
radius = (1.68/2)/12;    % diameter of 1.68 inches
mass = (1.62/16)/grav;   % mass of 1.62 ounces
rho = 0.0023769;         % density of air
area = pie*radius*radius;
inertia = 0.4*mass*radius*radius;   % inertia of a sphere

%  Launch conditions
v0 = [162 170 171 171 167 161 116 122 167 163]; % in mph
elev = [16.2 17.4 11.6 15.9 16.2 14.1 19.7 14.6 14.6 14.4]; % in deg
azim = [0.1 0.3 0.2 2.5 -0.1 5.4 -2.3 -0.5 2.1 -0.9];  % in deg, about Y (i.e. going left)
back = [5047 3349 4845 4762 3341 3340 3751 1684 2837 4184]; % in rpm, about global Z
side = [-2435 1063 421 618 -11 414 1146 655 1580 51];    % in rpm, about global Y
rifle = 0;   % in rpm, about global X

carry_sim = zeros([1 10]);
offline_sim = zeros([1 10]);

for i=1:10
    vx = (v0(i)*cos(elev(i)*pie/180)*cos(azim(i)*pie/180))*88/60;  % convert mph to ft/s
    vy = v0(i)*sin(elev(i)*pie/180)*88/60;
    vz = -(v0(i)*cos(elev(i)*pie/180)*sin(azim(i)*pie/180))*88/60;
    wx = rifle*pie/30;   % convert rpm to rad/s
    wy = side(i)*pie/30;
    wz = back(i)*pie/30;

    omega = sqrt(wx*wx + wy*wy + wz*wz);
    tx = wx/omega;
    ty = wy/omega;
    tz = wz/omega;

    t0 = 0;
    tf = 10;
    x0 = [vx, vy, vz, 0, 0, 0, omega, a, b, c, d, e]';   % launch conditions

    options = odeset('RelTol',1e-5,'AbsTol',1e-6);
    [t,x] = ode45('golf_eqns_2', [t0,tf], x0, options);

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

    carry = X(end);
    offline = Z(end);
    
    carry_sim(i) = carry;
    offline_sim(i) = offline;
end
 
    sums = zeros([1 10]);
    for i=1:10
        carry_diff = (carry_r(i) - carry_sim(i))^2;
        offline_diff = (offline_r(i) - offline_sim(i))^2;
        diff = sqrt(carry_diff + offline_diff);
        sums(i) = diff;
    end

    LPE_func = sum(sums)/10;
    
    LPE_func

end
