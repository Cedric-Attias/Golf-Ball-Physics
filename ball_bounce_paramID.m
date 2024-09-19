% Script that uses unconstrained and constrained NLP with least-squares to 
% identify model parameters for a model of continuous golf ball impacts
%
% Goal: Find the parameters for the Gonthier version of the Hunt-Crossley
% contact model for golf ball impacts (we will look at these types of
% models after reading week)
%

clc
clearvars
close all


% Let's make a fake data set to mock what we might see in doing our own
% experiments (can also see how sensitive our results are to noise)
cor_true = 0.84; % Let's say the true COR of our ball is constant (you'll check in Lab 1 as to whether this is true)
h1 = (0.5*1.47*1.47)/9.81; % Want our inbound speed to match Ismail and Stronge 2012
h2 = h1*cor_true*cor_true; % What the true bounce height

% Now we'll overlay some noise over random samples of the true values to
% simulate the act of data collection (for Lab 1 you don't need to sim data
% collection, we actually have data to use).
% - add variance assuming Gaussian additive  noise
num_samples = 10;
h1_variance = 0*0.0000001^2;
h2_variance = 0*0.0000001^2;
sampled_h1 = repmat(h1, [num_samples, 1]) + sqrt(h1_variance).*randn(num_samples,1);
sampled_h2 = repmat(h2, [num_samples, 1]) + sqrt(h2_variance).*randn(num_samples,1);

plot_exp_data = false;
if plot_exp_data
    subplot(1,2,1)
    hist(sampled_h1); %
    xlabel('Drop Height (m)');
    ylabel('Frequency')
    subplot(1,2,2)
    hist(sampled_h2);
    xlabel('Bounce Height (m)');
    ylabel('Frequency')
    input('');
    close all
    error;
end


%% Section A: Unconstrained parameter ID (least-squares)
% - use fmincon w/o specifying constraints to find the stiffness and damping coefficients (k and a) 
%   that match a set coefficient of restitution of a golfball, i.e. match
%   the drop experiment results

% First bit of code is just setting up necessary variables/functions
nlp_params.input_data = sampled_h1;
nlp_params.output_data = sampled_h2;
nlp_params.num_samples = num_samples;

model_params.m = 45/1000;
model_params.g = 9.81;
model_params.radius = 0.5*43/1000; 

% z0 = [50, 0.01];
z0 = [10e4,	1];
nlp_params.model_params = model_params;
cost_fun = @(z) cost(z, nlp_params);
options = optimoptions('fmincon', 'Algorithm','interior-point','Display','iter');


% Let's quick see how our results look with our initial guess (set
% check_initial_guess = true)
check_initial_guess = false;
if check_initial_guess
    model_params.k = z0(1);
    model_params.a = z0(2);
    
    ode_fun = @(t,X) ball_bounce_dynamics(t, X, model_params);
    event_fun = @(t,X) stop_ball_bounce(X, model_params);
    ode_options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Events', event_fun);
    
    tspan = [0,10];
    X0 = [h1, 0]; % assume no speed at t = 0
    [t_sim, X_sim] = ode45(ode_fun, tspan, X0, ode_options);
    animate_ball(t_sim, X_sim, model_params)
    disp(['Simulated bounce height: ', num2str(X_sim(end,1)), ' m']);
    input('');
    close all
    error();
end


% Now let's run the optimization and see how the results improve (set
% check_new_results to true)
check_new_results = false;
if check_new_results
   
    [opt_z, fval] = fmincon(cost_fun, z0, [], [], [], [], [], [], [], options);
    model_params.k = opt_z(1);
    model_params.a = opt_z(2);
    
    ode_fun = @(t,X) ball_bounce_dynamics(t, X, model_params); % Redefine each time since we're changing the model parameters
    event_fun = @(t,X) stop_ball_bounce(X, model_params);
    ode_options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Events', event_fun);
    
    tspan = [0,5];
    X0 = [h1, 0]; % assume no speed at t = 0
    [t_sim, X_sim] = ode45(ode_fun, tspan, X0, ode_options); 
    close all;
    animate_ball(t_sim, X_sim, model_params)
    disp(['Simulated bounce height: ', num2str(X_sim(end,1)), ' m']);
    input('');
    close all;

    % Does the optimization find a reasonable result? Let's check the max force
    % applied to the ball and duration the bodies are impacting each other
    [F_contact, t_contact] = calculation_contact_force(t_sim, X_sim, model_params);
    contact_duration = max(t_contact) - min(t_contact);
    plot(t_contact,F_contact);
    xlabel('Time (s)');
    ylabel('Normal Force (N)');
    disp(['Max Normal Force: ', num2str(max(F_contact)), ' N']);
    disp(['Contact Duration: ', num2str(contact_duration), ' s']);
    input('');
    close all;
    error;
end



% Ismail and Stronge 2012 max Force = 223.5 N; contact duration ~ 1 ms;
% See if we can get a better match by adding constraints to the problem
% --- let's try restricting the max contact force to the experimentally
% measured value using equality constraints
% --- we'll also use lb and ub to impose box inequality constraints (upper and lower
% limits) on the design variables
nlp_params.exp_normal_force = 223.5;
nonlcon_fun = @(z) nonlcon(z, nlp_params); 
lb = [0, 0];
ub = [10e9, 10e4];
check_constrained_results = false; % set to true to run constrained opt
if check_constrained_results

    [opt_z, fval] = fmincon(cost_fun, z0, [], [], [], [], lb, ub, nonlcon_fun, options);
    model_params.k = opt_z(1);
    model_params.a = opt_z(2);
    
    ode_fun = @(t,X) ball_bounce_dynamics(t, X, model_params); % Redefine each time since we're changing the model parameters
    event_fun = @(t,X) stop_ball_bounce(X, model_params);
    ode_options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Events', event_fun);
    
    tspan = [0,5];
    X0 = [h1, 0]; % assume no speed at t = 0
    [t_sim, X_sim] = ode45(ode_fun, tspan, X0, ode_options); 
    close all;
    animate_ball(t_sim, X_sim, model_params)
    disp(['Simulated bounce height: ', num2str(X_sim(end,1)), ' m']);
    input('');
    close all;

    % Does the optimization find a reasonable result? Let's check the max force
    % applied to the ball and duration the bodies are impacting each other
    [F_contact, t_contact] = calculation_contact_force(t_sim, X_sim, model_params);
    contact_duration = max(t_contact) - min(t_contact);
    plot(t_contact,F_contact);
    xlabel('Time (s)');
    ylabel('Normal Force (N)');
    disp(['Max Normal Force: ', num2str(max(F_contact)), ' N']);
    disp(['Contact Duration: ', num2str(contact_duration), ' s']);
    input('');
    close all;
    error;
    
end




% NB: Even with our final constrained problem, our Hunt-Crossley model is limited because in
% order to accurately capture the COR for different inbound conditions,
% need to adjust the "a" parameter relative to k and inbound velocity!!


%% Cost and Constraint Functions
function J = cost(z, nlp_params)

    % Step 1: Setup new model parameters by subbing new candidates for a
    % and k (from z vector) into our struct
    model_params = nlp_params.model_params;
    model_params.k = z(1);
    model_params.a = z(2);
    
    % Step 2: Simulation system using experimental inputs (in this case, our
    % input to the model is the drop/initial height)
    ode_fun = @(t,X) ball_bounce_dynamics(t, X, model_params);
    event_fun = @(t,X) stop_ball_bounce(X, model_params);
    ode_options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Events', event_fun);
    
    tspan = [0,5];
    num_samples = nlp_params.num_samples;
    input_data = nlp_params.input_data;
    sim_bounce_height = zeros(num_samples,1); % preallocate for speed
    for i = 1:num_samples
        X0 = [input_data(i,1), 0]; % assume no speed at t = 0
        [~, X_sim] = ode45(ode_fun, tspan, X0, ode_options);
        sim_bounce_height(i,1) = X_sim(end,1);
    end
    
    output_data = nlp_params.output_data;
    J = 1e2.*sum((sim_bounce_height-output_data).^2); % vectorized 

end



function [c,ceq] = nonlcon(z, nlp_params)
% Outputs: c = array of inequality constraints of form c(i) < 0 
%          ceq = array of equality constraints of form ceq(i) = 0


% Step 1: Setup new model parameters by subbing new candidates for a
% and k (from z vector) into our struct
model_params = nlp_params.model_params;
model_params.k = z(1);
model_params.a = z(2);

% Step 2: Simulation system using experimental inputs (in this case, our
% input to the model is the drop/initial height)
ode_fun = @(t,X) ball_bounce_dynamics(t, X, model_params);
event_fun = @(t,X) stop_ball_bounce(X, model_params);
ode_options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Events', event_fun);

tspan = [0,5];
num_samples = nlp_params.num_samples;
input_data = nlp_params.input_data;
peak_normal_force = zeros(num_samples,1); % preallocate for speed
for i = 1:num_samples
    X0 = [input_data(i,1), 0]; % assume no speed at t = 0
    [t_sim, X_sim] = ode45(ode_fun, tspan, X0, ode_options);
    [F_contact, ~] = calculation_contact_force(t_sim, X_sim, model_params);
    peak_normal_force(i) = max(F_contact);
end

% Try equality constraint?
% c = [];
% ceq = mean(peak_normal_force) - nlp_params.exp_normal_force;

% Try softer approach using inequality constraints
% c = 200 - mean(peak_normal_force);
% ceq = [];

% % Delete all
% c = [];
% ceq = [];

end




%% Extra Functions
function [F_contact, t_contact] = calculation_contact_force(t, X, model_params)
% returns array showing the force applied to the ball over the period of
% contact (F_contact) and an array of time over the period of contact
    k = model_params.k;
    a = model_params.a;
    radius = model_params.radius;
    
    T = length(t);
    F = zeros(T,1);
    for i = 1:T
        y = X(i,1);
        y_dot = X(i,2);
        delta = -(y-radius); 
        delta_dot = -y_dot;
        if delta > 0
            F(i) = k*delta*(3/2)*(1+a*delta_dot);
        else
            F(i) = 0;
        end
        
        if F(i) < 0 % enforce unilateral contact in the case of positive delta but large negative delta_dot
            F(i) = 0;
        end

    end

    contact_window = F > 0;
    t_contact = t(contact_window);
    F_contact = F(contact_window);

end


function animate_ball(t, X, model_params)
    
    fps = 60;
    t_resample = 0:(1/fps):t(end);
    [~, num_states] = size(X);
    X_resample = zeros(length(t_resample), num_states);
    for j = 1:num_states
        X_resample(:,j) = spline(t, X(:,j), t_resample);
    end
        
    for i = 1:length(t_resample)

        rect = findall(gcf,'Type', 'Rectangle'); 
        delete(rect); 
    
        y = X_resample(i,1);
        x = y.*0;
    
        r = model_params.radius;
        pos1 = [x-r y-r 2*r 2*r]; 
        rectangle('Position',pos1,'Curvature',[1 1], 'FaceColor', 'r')
        ylabel('Y (m)')
        xlabel('X (m)');

        axis([-0.15, 0.15, -0.10, 0.20]);
    %     xlim([x-0.1,x+2])
    %     ylim([0, 0.05]);
        set(gcf,'Position',[0 100 500 500])

        pause(1*(1/fps));
        drawnow;
    
    end

end



function [position,isterminal,direction] = stop_ball_bounce(X, params)
% Use subfunction here to capture event of ball landing during integration
  position = X(2); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = -1;  % -1 means that the X(2) at previous time_step should be > than the one at current iteration
end
    