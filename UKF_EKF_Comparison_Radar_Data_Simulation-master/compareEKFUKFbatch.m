clear all
close all

ProcessNoise_SigmaSquare = 1e-4;
MeasurementNoise_Sigma = [1e-3; 2e-2;3e-4];  % Range; Bearings
Correlated_Noise=8;

DataGenerate_GenerateTrajectory

store_ekf_errors = zeros(100, size(Observations, 2));
store_ukf_errors = zeros(100, size(Observations, 2));
store_ckf_errors = zeros(100, size(Observations, 2));
for iter = 1:100
% Data Generation
DataGenerate_GenerateTrajectory
DataGenerate_AddingMeasurementNoise

% Initialisation
[Q, R,D ] = EKFinit( ProcessNoise_SigmaSquare, MeasurementNoise_Sigma, T );
 Mean_ekf =[1980; 2020; 0; 0];
 Cov_ekf = diag([10; 10; 4; 4]);
store_ekf_errors(iter, 1) = norm(Mean_ekf(1:2, 1)-store_Target_real_state(1:2, 1));
[ Mean_ukf, Cov_ukf ] = UKFinit( [2010; 1990; 0; 0], diag([50; 50; 5; 5]), ProcessNoise_SigmaSquare, MeasurementNoise_Sigma, T,Correlated_Noise);
store_ukf_errors(iter, 1) = norm(Mean_ukf(1:2, 1)-store_Target_real_state(1:2, 1));
[ Mean_ckf, Cov_ckf ] = CKFinit( [2010; 1990; 0; 0], diag([50; 50; 5; 5]), ProcessNoise_SigmaSquare, MeasurementNoise_Sigma, T);
store_ckf_errors(iter, 1) = norm(Mean_ckf(1:2, 1)-store_Target_real_state(1:2, 1));
% Data Processing
for i = 2:size(Observations, 2)
%     [ Mean_ekf, Cov_ekf ] = EKFrun( Mean_ekf, Cov_ekf, Observations(:, i), Q, R, T );
%     store_ekf_errors(iter, i) = norm(Mean_ekf(1:2, 1)-store_Target_real_state(1:2, i));
    [ Mean_ukf, Cov_ukf ] = UKFrun_myself( Mean_ukf, Cov_ukf, Observations(:, i), Q, R, T);
    store_ukf_errors(iter, i) = norm(Mean_ukf(1:2, 1)-store_Target_real_state(1:2, i));
    [ Mean_ckf, Cov_ckf ] = CKFrun( Mean_ckf, Cov_ckf, Observations(:, i), Q, R, T);
    store_ckf_errors(iter, i) = norm(Mean_ckf(1:2, 1)-store_Target_real_state(1:2, i));
end
fprintf('Monte Carlo Simulation: %d Finished!... \n', iter);
end

ekf_final = mean(store_ekf_errors, 1);
ukf_final = mean(store_ukf_errors, 1);
ckf_final = mean(store_ckf_errors, 1);
figure(1);
hold on;
plot(ekf_final, 'r-');
plot(ukf_final, 'b-');
plot(ckf_final, 'g-');
legend('ekf', 'ukf','ckf');
xlabel('Timestep');
ylabel('Euclidean Error');    