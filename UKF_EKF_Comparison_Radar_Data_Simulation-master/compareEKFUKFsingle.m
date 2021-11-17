close all
clear
ProcessNoise_SigmaSquare = 1e-4;
MeasurementNoise_Sigma = [1e-3; 2e-2; 3e-4];  % Range; Bearings
Correlated_Noise=8;

% Plot
figure(1);
hold on;

% Data Generation
DataGenerate_GenerateTrajectory

% Initialisation
store_ekf_errors = zeros(1, size(Observations, 2));
store_ukf_errors = zeros(1, size(Observations, 2));
store_ckf_errors = zeros(1, size(Observations, 2));
store_sckf_errors = zeros(1, size(Observations, 2));
[Qe, Re,De] = EKFinit(  ProcessNoise_SigmaSquare, MeasurementNoise_Sigma, Correlated_Noise  );
Mean_ekf =[1980; 2020; 0; 0];
Cov_ekf = diag([10; 10; 4; 4]);
store_ekf_errors(1, 1) = norm(Mean_ekf(1:2, 1)-store_Target_real_state(1:2, 1));
[ Mean_ukf, Cov_ukf,Qu,Ru ,Du] = UKFinit( [1980; 2020; 0; 0], diag([10; 10; 4; 4]), ProcessNoise_SigmaSquare, MeasurementNoise_Sigma ,T,Correlated_Noise);
store_ukf_errors(1, 1) = norm(Mean_ukf(1:2, 1)-store_Target_real_state(1:2, 1));
[ Mean_ckf, Cov_ckf,Qc,Rc] = CKFinit( [1980; 2020; 0; 0], diag([10; 10; 4; 4]), ProcessNoise_SigmaSquare, MeasurementNoise_Sigma ,T);
store_ckf_errors(1, 1) = norm(Mean_ckf(1:2, 1)-store_Target_real_state(1:2, 1));
[ Mean_sckf, Cov_sckf,Qsc,Rsc,S_input] = SCKFinit( [1980; 2020; 0; 0], diag([10; 10; 4; 4]), ProcessNoise_SigmaSquare, MeasurementNoise_Sigma ,T);
store_ckf_errors(1, 1) = norm(Mean_ckf(1:2, 1)-store_Target_real_state(1:2, 1));
plot(Mean_ekf(1), Mean_ekf(2), 'r.');
plot(Mean_ukf(1), Mean_ukf(2), 'b.');
plot(Mean_ckf(1), Mean_ckf(2), 'g.');
plot(Mean_sckf(1), Mean_sckf(2), 'c.');
plot(store_Target_real_state(1, 1), store_Target_real_state(2, 1), 'g.');
% Data Processing
for i = 2:size(Observations, 2)
    [ Mean_ekf, Cov_ekf ] = EKFrun( Mean_ekf, Cov_ekf, Observations(:, i), Qe, Re, T );
     store_ekf_errors(1, i) = norm(Mean_ekf(1:2, 1)-store_Target_real_state(1:2, i));

    [ Mean_ukf, Cov_ukf ] = UKFrun_myself( Mean_ukf, Cov_ukf, Observations(:, i), Qu, Ru, T );
    store_ukf_errors(1, i) = norm(Mean_ukf(1:2, 1)-store_Target_real_state(1:2, i));
   
    [ Mean_ckf, Cov_ckf ] = CKFrun( Mean_ckf, Cov_ckf, Observations(:, i), Qc, Rc, T );
    store_ckf_errors(1, i) = norm(Mean_ckf(1:2, 1)-store_Target_real_state(1:2, i));   

    plot(Mean_ekf(1), Mean_ekf(2), 'r.');
    plot(Mean_ukf(1), Mean_ukf(2), 'b.');
    plot(Mean_ckf(1), Mean_ckf(2), 'g.');
    plot(Mean_sckf(1), Mean_sckf(2), 'c.');
    plot(store_Target_real_state(1, i), store_Target_real_state(2, i), 'y.');
end
figure(2);
hold on;
xlabel('Timestep');
ylabel('Euclidean Error');
plot(store_ekf_errors, 'r-');
plot(store_ukf_errors, 'b-');
plot(store_ckf_errors, 'g-');
plot(store_sckf_errors, 'c-');
legend('EKF', 'UKF','CKF');
