clear;
clc;
close all;

%% Part 1: Experimental data and stiffness search range

% Experimental gauge pressure, unit: kPa  
% (= Measured pressure - atmosphere pressure)
pressure_exp_kPa = [0; 2; 3; 4];

% Experimental angle, unit: degree
theta_exp_deg = [46; 161; 178; 180];

% Make sure both are column vectors
pressure_exp_kPa = pressure_exp_kPa(:);
theta_exp_deg = theta_exp_deg(:);

% Check the number of experimental points
if length(pressure_exp_kPa) ~= length(theta_exp_deg)
    error('Pressure data and angle data must have the same length.');
end

% Rotational spring stiffness search range
K_min = 0.04;
K_max = 0.09;

% Number of stiffness candidates
num_K = 10;

% Generate stiffness candidates
K_candidates = linspace(K_min, K_max, num_K);

%% Part 2: Run simulations

num_pressure = length(pressure_exp_kPa);

% Each column corresponds to one stiffness
% Each row corresponds to one experimental pressure
theta_sim_all = nan(num_pressure,num_K);

% Record whether each stiffness completed successfully
simulation_success = false(num_K,1);

for i = 1:num_K

    K_current = K_candidates(i);

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('Running stiffness %d of %d\n',i,num_K);
    fprintf('Current K = %.8f\n',K_current);
    fprintf('========================================\n');

    current_success = true;

    for j = 1:num_pressure

        pressure_current = pressure_exp_kPa(j);

        fprintf('  Pressure %d of %d: %.4f kPa\n', ...
            j,num_pressure,pressure_current);

        try
            theta_sim_all(j,i) = ...
                run_structure_calibration( ...
                K_current, ...
                pressure_current);

            if ~isfinite(theta_sim_all(j,i))
                warning( ...
                    'K = %.8f, pressure = %.4f kPa returned NaN or Inf.', ...
                    K_current,pressure_current);

                current_success = false;
                break
            end

        catch ME
            warning( ...
                'Simulation failed for K = %.8f at pressure %.4f kPa.\n%s', ...
                K_current,pressure_current,ME.message);

            current_success = false;
            break
        end

    end

    simulation_success(i) = current_success;

end

%% Part 3: Calculate RMSE

RMSE_values = inf(num_K,1);

for i = 1:num_K

    if simulation_success(i)

        % residual = theta_sim_all(:,i)-theta_exp_deg;
        % 
        % RMSE_values(i) = sqrt(mean(residual.^2));

        residual = theta_sim_all(:,i)-theta_exp_deg;

        base_RMSE = sqrt(mean(residual.^2));
        
        % Constraint on the second experimental point:
        % simulation angle at 2 kPa cannot exceed experiment
        excess_second_point = max(0,theta_sim_all(2,i)-theta_exp_deg(2));
        
        penalty_weight = 10;
        
        RMSE_values(i) = base_RMSE + ...
        penalty_weight*excess_second_point;

    else

        RMSE_values(i) = Inf;

    end

end

[best_RMSE,best_index] = min(RMSE_values);

if ~isfinite(best_RMSE)
    error('All stiffness candidates failed. No valid best stiffness was found.');
end

best_K = K_candidates(best_index);

best_theta_sim_deg = theta_sim_all(:,best_index);

fprintf('\n');
fprintf('========================================\n');
fprintf('Calibration completed\n');
fprintf('Best stiffness = %.10f\n',best_K);
fprintf('Minimum RMSE = %.6f degree\n',best_RMSE);
fprintf('========================================\n');

calibration_result = table( ...
    K_candidates(:), ...
    RMSE_values, ...
    simulation_success, ...
    'VariableNames', ...
    {'Stiffness','RMSE_degree','SimulationSuccess'});

disp(calibration_result)

%% Part 4: Plot best simulation versus experiment

figure('Color','w');

plot( ...
    pressure_exp_kPa, ...
    theta_exp_deg, ...
    'o-', ...
    'LineWidth',1.8, ...
    'MarkerSize',8);

hold on;

plot( ...
    pressure_exp_kPa, ...
    best_theta_sim_deg, ...
    's--', ...
    'LineWidth',1.8, ...
    'MarkerSize',8);

grid on;
box on;

xlabel('Gauge pressure (kPa)');
ylabel('Angle (degree)');

title(sprintf( ...
    'Best fit: K = %.6g, RMSE = %.3f degree', ...
    best_K,best_RMSE));

legend( ...
    'Experiment', ...
    'Simulation', ...
    'Location','best');

figure('Color','w');

plot( ...
    K_candidates, ...
    RMSE_values, ...
    'o-', ...
    'LineWidth',1.5, ...
    'MarkerSize',7);

hold on;

plot( ...
    best_K, ...
    best_RMSE, ...
    'p', ...
    'MarkerSize',14, ...
    'LineWidth',2);

grid on;
box on;

xlabel('Rotational spring stiffness');
ylabel('RMSE (degree)');

title('RMSE versus rotational spring stiffness');

legend( ...
    'Tested stiffness', ...
    'Best stiffness', ...
    'Location','best');