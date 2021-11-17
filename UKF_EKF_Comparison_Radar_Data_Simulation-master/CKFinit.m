function [ Mean, Cov, Q, R] = CKFinit( InitX, InitCov, ProcessNoise_Sigma, MeasurementNoise_Sigma,T)
%   Detailed explanation goes here
Q = ProcessNoise_Sigma * ...
    [(1/3)*(T^3), 0,            (1/2)*(T^2), 0;
     0,           (1/3)*(T^3),  0,           (1/2)*(T^2);
     (1/2)*(T^2),  0,           T,            0;
     0,            (1/2)*(T^2), 0,            T];
R = diag(MeasurementNoise_Sigma) .^ 2;
Mean = InitX;
Cov = InitCov;
end

%(1/3)*(T^3) = 1/3*T^3
%  ProcessNoise_SigmaSquare = 1e-4;
%  MeasurementNoise_Sigma = [1e-3; 2e-2]; 