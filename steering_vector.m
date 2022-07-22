function [a] = steering_vector(M, theta, Delta_f, mu_0, fc, d, Ts)      % M = number of sensor, L = signal length                                                              % theta = DOA (degree), freq = frequency of single tone signal
    
    c = 3 * 10^8;
    % Propagation delay: tau = D + kappa; % D: integer delay % kappa: fractional delay
    tau = (0:M-1).' * d * sin(theta)/c; %  a (M, 1) vector
    D = round(tau/Ts) * Ts;
    kappa = tau - D;

    %% Generate $a(\theta, \Delta_f)$
    a = exp(-1j * 2 * pi * fc * tau) .* exp(-1j * 2 * pi * Delta_f * kappa) .* exp(1j * pi * mu_0 * kappa .^ 2);
end
