function [ Mean_Update, Cov_Update ] = UKFrun_ZengGuang( Mean, Cov, Observation, Q, R, T )
%position and velocity transfer matrix
F = [1, 0, T, 0;
     0, 1, 0, T;
     0, 0, 1, 0;
     0, 0, 0, 1];
% Parameters for UKF
n = size(Mean, 1);
Alpha = 0.25;
K = 1;
Beta = 4;
Lambda = Alpha^2*(n+K)-n;
C = sqrt(n);
SPs = Find_SigmaPoints(Mean, Cov, C );%sigma points
NumOfSP = size(SPs, 2);%number of sigma points
%weights
W_Mean = [Lambda/(2*n+Lambda), repmat(1/(2*(2*n+Lambda)), [1, 2*n])];% 1*(2n+1) 
W_Cov = [Lambda/(2*n+Lambda)+(1-Alpha^2+Beta), repmat(1/(2*(n+Lambda)), [1, 2*n])];

Mean_predict_star=F*SPs;
Weighted_x = Mean_predict_star * diag(W_Mean); %Mean_predict_star  : n*(2n+1)  /  diag(W_mean)  : (2n+1) * (2n+1)
Mean_predict = sum(Weighted_x, 2);
Mu_x_subtract = Mean_predict_star - repmat(Mean_predict, [1, NumOfSP]);
Cov_predict= Mu_x_subtract * diag(W_Cov) * Mu_x_subtract'+Q;

SPs2 = Find_SigmaPoints(Mean_predict, Cov   _predict, C);

NumOfSP2 = size(SPs2, 2);

Mu_x_subtract2 = SPs2 - repmat(Mean_predict, [1, NumOfSP2]);
SigmaPoints_Obs = zeros(2, NumOfSP2);
for i = 1:NumOfSP2
    SigmaPoints_Obs(:, i) = h(SPs2(:, i));
end
Weighted_y = SigmaPoints_Obs * diag(W_Mean);
Mu_y = sum(Weighted_y, 2);
Mu_y_subtract = SigmaPoints_Obs - repmat(Mu_y, [1, NumOfSP2]);

P_yy= Mu_y_subtract * diag(W_Cov) * Mu_y_subtract'+R;
P_xy = Mu_x_subtract2 * diag(W_Cov) * Mu_y_subtract';
K = P_xy*inv(P_yy);
Mean_Update = Mean_predict + K*(Observation - Mu_y);%Mu_x: Mean u_x
Cov_Update = Cov_predict-P_xy*inv(P_yy)*P_xy';

end

function [XinPolar] = h(XinCartesian)
    x = XinCartesian(1); %第一行那一个数 position x
    y = XinCartesian(2);  %第二行那个数  position y
    XinPolar = [sqrt(x^2+y^2); atan2(y, x)]; % 距离 还有 方向角
end

function [SPs] = Find_SigmaPoints(XinCartesian, Cov, C)
%Using eigenvalue decomposition 特征值分解，奇异值分解（SVD）
NumOfDim = size(XinCartesian, 1);%维度是n
NumOfPoints = NumOfDim*2;%2n
[V, D] = eig(Cov);

pnts = zeros(NumOfDim, NumOfPoints);%n行  2n列 pnts  points
for i = 1:NumOfDim
    pnts(i, i) = sqrt(abs(D(i, i)));
    pnts(i, i+NumOfDim) = -sqrt(abs(D(i, i)));
end
% scale samples such that they have the same second order states as a standard normal
Covariance_pnts = C*pnts;
Rotated_pnts = V*Covariance_pnts;
SigmaPoints = repmat(XinCartesian, [1, NumOfPoints]) + Rotated_pnts;%n*2n  行列数
SPs = [XinCartesian, SigmaPoints]; %所以并没有舍弃初始值， 最后的数目确实是2n+1个
end