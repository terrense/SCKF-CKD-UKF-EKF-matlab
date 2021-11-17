function [ Mean_Update, Cov_Update ] = EKFrun( Mean, Cov, Observation, Q, R, T )
%EFKSIMULATE Summary of this function goes here
%   Detailed explanation goes here
F = [1, 0, T, 0;
     0, 1, 0, T;
     0, 0, 1, 0;
     0, 0, 0, 1];
% Predict
Mean_Predict = F * Mean;
Cov_Predict = F*Cov*F' + Q;
% Update
H = Calculate_Jacobian(Mean_Predict);
S = H * Cov_Predict * H' + R;
% [U,D] = modchol(S);
% X1=	H'*inv(U');
% X2=X1*inv(D);
% X=X2*inv(U);
% K=Cov_Predict*X;
K = Cov_Predict * H' * inv(S);
Mean_Update = Mean_Predict + K* (Observation - h(Mean_Predict));
Cov_Update = (eye(4) - K*H) * Cov_Predict;
end
function [XinPolar] = h(XinCartesian)
  x = XinCartesian(1); %��һ����һ���� position x
  y = XinCartesian(2);  %�ڶ����Ǹ���  position y
  Vx= XinCartesian(3);
  Vy= XinCartesian(4);
XinPolar = [sqrt(x^2+y^2); atan2(y, x);sqrt(Vx^2+Vy^2)]; % 
end
function [J] = Calculate_Jacobian(XinCartesian)
    x = XinCartesian(1);
    y = XinCartesian(2);
    Vx= XinCartesian(3);
    Vy= XinCartesian(4);
    J = [x/sqrt(x^2+y^2),    y/sqrt(x^2+y^2),   0,                  0;
        -y/(x^2+y^2),       x/(x^2+y^2),        0,                  0;
         0,                  0,                 Vx/sqrt(Vx^2+Vy^2), Vy/sqrt(Vx^2+Vy^2);];
end