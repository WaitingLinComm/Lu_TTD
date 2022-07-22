%% Workspace Initialization.
clc; clear; close all;
%% Parameters.
fc   = 1 * 10^9;
fmax = 1.04 * 10^9;     % (GHz)
fmin = 0.96 * 10^9;     % (GHz)
B_F  = 0.08;            % Fractional Bandwidth(FB)
B    = B_F * fc;        % B = 0.1(GHz)
Tp   = 40 * 10^(-6);    % 40($\mu$ s) : Pulse duration
fs   = 10 * B;          % Sampling frequency
Ts   = 1/fs;            % Sampling period
c    = 3 * 10^8;
theta_0     = 20;                 % (deg) DOA of desired signal
theta_0_rad = (theta_0/180) * pi; % (rad) DOA of desired signal
theta_q     = (-65:0.5:-55);      % (deg) 想壓抑的 sidelobe 角度區間
theta_q_rad = (theta_q/180) * pi; % (rad) 想壓抑的 sidelobe 角度區間

Delta_f = [-B/2 0 B/2]; % Delta_f = {-B/2, 0 , B/2} % Steering vector中的參數
mu_0 = B/Tp;            % chirp rate: parameters of LFM signal
%% ULA
M          = 16;           % Consider a ULA with M isotropic antenna
lambda_min = c/fmax;
d          = lambda_min/2;    
%d = lambda_min/(2 * 1.15); % 此 paper 非使用 d = lambda_min/2 = c/(2 * fmax)   
%% Initializaiton
rho_1    = {};
rho_1{1} = 1e-4; %rho_1{1} = 1e-10; % p.3446
rho_2    = 0.01; % not mentioned in paper
xi_1     = 0.999;
xi_2     = 1.001;
varsigma = 1e-4;

a_theta0_Deltaf0 = steering_vector(M, theta_0_rad, 0, mu_0, fc, d, Ts);
w_0              = a_theta0_Deltaf0/sqrt(M);
w                = {};
w{1}             = w_0;

u            = {};
u{1}         = 1e-3; % not mentioned in paper %u{1} = rand(1, 1);
eta          = {};
eta_prime    = {};
eta_prime{1} = 0.001;
Delta_r      = {};
Delta_r{1}   = 0.001; % initial Delta_r should >= varsigma

%Gamma = 1; % minimum mainlobe level constraint
Gamma = 10 * log10(20 * log10(a_theta0_Deltaf0' * a_theta0_Deltaf0)/sqrt(M) - 0.6); % P = 2.7093 %From paper
%% CA algorithm
Kmax = 30; % Kmax = 2 * 10^4
for k = 1:Kmax
    display(['Iteration: ', num2str(k)])
    if Delta_r{end} <= varsigma
        break
    else
            %% Step 1: update rho_1
        if k >= 2
            if u{end} <= xi_1 * u{end - 1}
                rho_1{end + 1} =  rho_1{end};
            else
                rho_1{end + 1} =  xi_2 * rho_1{end};
            end
        end
            %% Step 2
        cvx_begin quiet
        variable eta_k; % eta_k is a real value 
        variable u_k;   % u_k is a real value 
        variable w_k(M, 1) complex;

        objective_function = eta_k + rho_1{end} * u_k + rho_2 * (w_k - w{end})' * (w_k - w{end}); % sol1
        %objective_function = eta_k + rho_1{end} * u_k + rho_2 * power(2, norm(w_k - w{end}, 2)); % sol2

        minimize objective_function
        subject to 
            for q = 1:length(theta_q_rad)
                for p = 1:length(Delta_f)
                    a_thetaq_Deltafp = steering_vector(M, theta_q_rad(q), Delta_f(p), mu_0, fc, d, Ts);
                    R_qp = a_thetaq_Deltafp * a_thetaq_Deltafp';
                    %a_thetaq_Deltafp_normalized = a_thetaq_Deltafp/norm(a_thetaq_Deltafp, 2); % Normalized steering vector
                    %R_qp = a_thetaq_Deltafp_normalized * a_thetaq_Deltafp_normalized';        % with Normalized steering vector                   
                    w_k' * R_qp *  w_k <= eta_k;
                end   
            end
            real(w_k' * a_theta0_Deltaf0) >= sqrt(Gamma);
            imag(w_k' * a_theta0_Deltaf0) == 0;
            w_k' * w_k <= 1;
            1 - u_k - w{end}' * w{end} - real(2 * w{end}' * (w_k - w{end})) <= 0;
            u_k >= 0;
            u_k <= u{end};
        cvx_end
        u{end + 1} = u_k;
        w{end + 1} = w_k;
        eta{end + 1} = eta_k;
            %% Step 3: Stop criterion
        eta_prime{end + 1} = eta{end} + rho_1{end} * u{end};
        Delta_r{end + 1} = abs((eta_prime{end} - eta_prime{end - 1})/eta_prime{end - 1});
    end % end if Delta_r <= varsigma

    w{end + 1} = w_k; % check: w_k' * w_k = 1
    eta{end + 1} = eta_k;
end
w_opt = w{end};

%% Save file
%save('w_opt.mat','w_opt');
%% Plot Beampattern. (Please refer to main_plot_beampattern.m)
theta_grid = (-90 : 0.1 : 90);          % Plotted theta_grid interval(in degrees).
theta_grid_rad = theta_grid./180 .* pi; %rad
Ntheta = length(theta_grid);            % Ntheta = 1801
%f_range = [0];
%f_range = [-10:10];
f_range = [-10^3:10^3];
 
for q = 1:length(theta_grid_rad)
    for p = 1:length(f_range)
        a_thetaq_Deltafp = steering_vector(M, theta_grid_rad(q), f_range(p), mu_0, fc, d, Ts);
        a_thetaq_Deltafp_normalized = a_thetaq_Deltafp/norm(a_thetaq_Deltafp, 2);
        %Beam_Pat(p, q) = w_0' * a_thetaq_Deltafp;
        %Beam_Pat(p, q) = w_opt' * a_thetaq_Deltafp;
        Beam_Pat(p, q) = w_opt' * a_thetaq_Deltafp_normalized;
    end   
end
%{
figure
plot(theta_grid, (abs(Beam_Pat).^2));

figure
plot(theta_grid, 10 * log10(abs(Beam_Pat).^2));
title('Beampattern');
ylabel('Magnitude response(dB)');
xlabel('$Angle: \theta$ (deg)','interpreter','latex');
xlim([-90,90]);
%}


figure
mesh(theta_grid, f_range, 10 * log10(abs(Beam_Pat).^2));
title('Beampattern');
ylabel('Frequency: $f$','interpreter','latex');
xlabel('$DOA: \theta$ (deg)','interpreter','latex');
zlabel('$|B(f, \theta)|(dB)$','interpreter','latex');

