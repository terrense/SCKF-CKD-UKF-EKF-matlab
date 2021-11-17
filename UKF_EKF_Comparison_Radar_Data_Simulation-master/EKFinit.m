function [Q, R,D] = EKFinit(ProcessNoise_Sigma, MeasurementNoise_Sigma, Correlated_Noise )
%   Detailed explanation goes here
T=1;
Q = ProcessNoise_Sigma * ...
    [1/3*T^3, 0,       1/2*T^2, 0;
     0,       1/3*T^3, 0,       1/2*T^2;
     1/2*T^2, 0,       T,       0;
     0,       1/2*T^2, 0,       T];
R = diag(MeasurementNoise_Sigma) .^ 2;
% Mean = InitX;
% Cov = InitCov;
D=[Correlated_Noise ,Correlated_Noise 0,0;
    0,0,0,0]';
end

