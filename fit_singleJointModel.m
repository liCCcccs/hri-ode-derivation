


load hri_data.mat


% Sample data (replace with your actual data)
th = hri_data.ang_diff;
dth = hri_data.vel_diff;
ang = rad2deg(hri_data.r_q5);
tau_int = hri_data.tau_h_eqv;

% Create the design matrix A
A = [th, dth, ang];

% Perform a linear regression to fit the model k = ax + by + cz
coefficients = A\tau_int;

% Extract the coefficients
a = coefficients(1);
b = coefficients(2);
c = coefficients(3);

tau_int_est = a * th + b * dth + c * ang;

% Display the coefficients
fprintf('a = %.4f\n', a);
fprintf('b = %.4f\n', b);
fprintf('c = %.4f\n', c);

figure('Renderer', 'painters', 'Position', [300 300 800 800])
hold on
plot(th, tau_int)
plot(th, tau_int_est)
legend("tau actual", "tau estimate")

r_squared = get_r2(tau_int, tau_int_est)


function [r_squared] = get_r2(actual, predicted)
    % Calculate the mean of actual values
    mean_actual = mean(actual);
    
    % Calculate the total sum of squares (TSS)
    tss = sum((actual - mean_actual).^2);
    
    % Calculate the residual sum of squares (RSS)
    rss = sum((actual - predicted).^2);
    
    % Calculate R-squared (R2)
    r_squared = 1 - (rss / tss);
end







