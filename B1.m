clear;
clc;
N_values = [1000, 2000];
%Output table with calculated optimal results. 
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('|   N   |  Mean[y]  |  Var[y]  | AR(1) a1 | Var[E1] | AR(2) a1 | AR(2) a2 | Var[E2] |\n');
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('|  inf  |  0.0000   |  4.3330  | -0.2690  | 4.0200  | -0.2512  |  0.0670  |  4.0012 |\n');
    fprintf('-------------------------------------------------------------------------------------\n');
    fprintf('1 Set of data\n');
    fprintf('-------------------------------------------------------------------------------------\n');                       

    %recursive algorithm for data generation
for k = 1:length(N_values)
    N = N_values(k);
    y = zeros(1, N);
    e = randn(1, N) * 2; 
    for t = 2:N
        y(t) = (-1/2) * y(t-1) + e(t) + (1/4) * e(t-1);
    end
    
    y = y(:);
    mean_y_1set = mean(y); %expected value of stationary process
    variance_y_1set = var(y); %var of stationary process

    %AR(1)
    Ar1_y_lag1 = [0; y(1:end-1)]; % y(t-1)
    theta1_est_1set = sum(Ar1_y_lag1 .* y) / sum(Ar1_y_lag1 .* Ar1_y_lag1);
    
    error1 = y - Ar1_y_lag1 * theta1_est_1set;%error
    variance_of_error1 = var(error1);%variance of the error
    
    %AR(2)
    Ar2_y_lag1 = [0; y(1:end-1)]; % y(t-1) 
    Ar2_y_lag2 = [0; 0; y(1:end-2)]; % y(t-2)
    theta2_est_1set = [sum(Ar2_y_lag1 .* Ar2_y_lag1), sum(Ar2_y_lag1 .* Ar2_y_lag2); ...
        sum(Ar2_y_lag1 .* Ar2_y_lag2), sum(Ar2_y_lag2 .* Ar2_y_lag2)] \ [sum(Ar2_y_lag1 .* y); sum(Ar2_y_lag2 .* y)];
    % Calculate the error
    error2 = y - Ar2_y_lag1 *theta2_est_1set(1) - Ar2_y_lag2*theta2_est_1set(2);
    % Calculate the variance of the error
    variance_of_error2 = var(error2);

    fprintf('|%6s | %7.4f   | %7.4f  | %7.4f  |%7.4f  | %7.4f  | %7.4f  | %7.4f |\n', num2str(N),...
        mean_y_1set, variance_y_1set, theta1_est_1set, var(error1), theta2_est_1set(1), theta2_est_1set(2), var(error2));


end

fprintf('-------------------------------------------------------------------------------------\n');                        
fprintf('200 Sets of data\n');
fprintf('-------------------------------------------------------------------------------------\n');                        


N_values = [1000, 2000];%loop

% Initialize matrices for storing the values for every iteration of the loop.
theta1_estimates = zeros(200, length(N_values));
theta2_estimates = zeros(200, 2, length(N_values));
mean_y_estimates = zeros(200, length(N_values));
variance_y_estimates = zeros(200, length(N_values));
Var_Error200_1 = zeros(200, length(N_values));
Var_Error200_2 = zeros(200, length(N_values));

for k = 1:length(N_values)
    N = N_values(k);
    % Loop to generate 200 sets of data
    for i = 1:200
        
        y = zeros(1, N);
        e = randn(1, N) * 2;
        for t = 2:N
            y(t) = (-1/2) * y(t-1) + e(t) + (1/4) * e(t-1);
        end
        y = y(:);

        %AR(1)
        Ar1_y_lag1 = [0; y(1:end-1)]; % y(t-1)
        theta1_est_200set = sum(Ar1_y_lag1 .* y) / sum(Ar1_y_lag1 .* Ar1_y_lag1);
        Error200_1 = y - Ar1_y_lag1 * theta1_est_1set;

        %AR(2)
        Ar2_y_lag1 = [0; y(1:end-1)]; % y(t-1) 
        Ar2_y_lag2 = [0; 0; y(1:end-2)]; % y(t-2)
        theta2_est_200set = [sum(Ar2_y_lag1 .* Ar2_y_lag1), sum(Ar2_y_lag1 .* Ar2_y_lag2);...
            sum(Ar2_y_lag1 .* Ar2_y_lag2), sum(Ar2_y_lag2 .* Ar2_y_lag2)] \ [sum(Ar2_y_lag1 .* y); sum(Ar2_y_lag2 .* y)];
        Error200_2 = y - Ar2_y_lag1 *theta2_est_1set(1) - Ar2_y_lag2*theta2_est_1set(2);

        % Store the estimates in the corresponding matrices
        mean_y_estimates(i, k) = mean(y);
        variance_y_estimates(i, k) = var(y);
        Var_Error200_1(i, k)=var(Error200_1);
        Var_Error200_2(i, k)=var(Error200_2);
        theta1_estimates(i, k) = theta1_est_200set;
        theta2_estimates(i, :, k) = theta2_est_200set;
    end
    %add the results to the output table
    fprintf('|%6s | %7.4f   | %7.4f  | %7.4f  |%7.4f  | %7.4f  | %7.4f  | %7.4f |\n', ...
        num2str(N), mean(mean_y_estimates(:, k)), mean(variance_y_estimates(:, k)), mean(theta1_estimates(:, k)), ...
        Var_Error200_1(i, k), mean(theta2_estimates(:, 1, k)), mean(theta2_estimates(:, 2, k)), Var_Error200_2(i, k));
end
fprintf('-------------------------------------------------------------------------------------\n');                       


% Mat array to store the squared differences for each N
squared_diff = zeros(200, length(N_values));
% Mat to store the mean squared differences for each N
mean_squared_diff_theta1 = zeros(1, length(N_values));

for k = 1:2
    N = N_values(k);
    mean_theta1 = mean(theta1_estimates(:, k));%mean of the theta1 estimates
    % Calculate the squared differences for this N
    for i = 1:200
        squared_diff(i, k) = (theta1_estimates(i, k) - mean_theta1)^2;
    end
    mean_squared_diff_theta1(k) = mean(squared_diff(:, k));% mean of the squared differences for this N
end

% Initialize matrices to store the mean squared differences for theta2 parameters and their covariance
mean_squared_diff_theta2 = zeros(2, 2, length(N_values));
cov_theta2 = zeros(1, length(N_values));  % This will store the covariance between a1 and a2 of theta2


for k = 1:length(N_values)
    N = N_values(k);

    % Compute the means of the theta2 estimates for this N
    mean_theta2_a1 = mean(theta2_estimates(:, 1, k));
    mean_theta2_a2 = mean(theta2_estimates(:, 2, k));
    % Initialize arrays for squared differences and products for theta2
    squared_diff_theta2_a1 = zeros(200, 1);
    squared_diff_theta2_a2 = zeros(200, 1);
    product_diff_theta2 = zeros(200, 1);

    % Calculate the squared differences and product differences for theta2
    for i = 1:200
        squared_diff_theta2_a1(i) = (theta2_estimates(i, 1, k) - mean_theta2_a1)^2;
        squared_diff_theta2_a2(i) = (theta2_estimates(i, 2, k) - mean_theta2_a2)^2;
        product_diff_theta2(i) = (theta2_estimates(i, 1, k) - mean_theta2_a1) * (theta2_estimates(i, 2, k) - mean_theta2_a2);
    end

    % Compute the mean of the squared differences for theta2, which is the variance
    mean_squared_diff_theta2(1, 1, k) = mean(squared_diff_theta2_a1);
    mean_squared_diff_theta2(2, 2, k) = mean(squared_diff_theta2_a2);

    % Compute the mean of the product differences for theta2, which is the covariance
    cov_theta2(k) = mean(product_diff_theta2);
end

% Add the variance-covariance matrix in the output table
fprintf('AR(1) Empirical Variance:\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('Empirical Variance of Θ1 for N=1000:\n');
fprintf('| %10.4f |\n', mean_squared_diff_theta1(1));
fprintf('Empirical Variance of Θ1 for N=2000:\n');
fprintf('| %10.4f |\n', mean_squared_diff_theta1(2));
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('AR(2) Empirical Variance:\n');
fprintf('-------------------------------------------------------------------------------------\n');
for k = 1:length(N_values)
    var_theta2=[mean_squared_diff_theta2(1, 1, k) cov_theta2(k); cov_theta2(k) mean_squared_diff_theta2(2, 2, k)];
    fprintf('Empirical Variance of Θ2 for N=%d:\n', N_values(k));
    fprintf('| %10.4f %10.4f |\n', var_theta2(1,1), var_theta2(1,2));
    fprintf('| %10.4f %10.4f |\n', var_theta2(2,1), var_theta2(2,2));
end
