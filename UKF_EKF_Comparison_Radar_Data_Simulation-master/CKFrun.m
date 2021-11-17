function [ Mean_Update, Cov_Update ] = CKFrun( Mean, Cov, Observation, Q, R, T )
%position and velocity transfer matrix
F = [1, 0, T, 0;
     0, 1, 0, T;
     0, 0, 1, 0;
     0, 0, 0, 1];
% Parameters for CKF
n = size(Mean, 1);
%weights
Weight=repmat(1/(2*n),[1,2*n]);
%Time Update
SPs1 = Find_SigmaPoints(Mean, Cov);%sigma points     (n X 2n)
NumOfSP = size(SPs1, 2);%number of sigma points
Mean_predict_star=F*SPs1;
Weighted_x = Mean_predict_star * diag(Weight);
Mean_predict = sum(Weighted_x, 2);
Mu_x_subtract = Mean_predict_star - repmat(Mean_predict, [1, NumOfSP]);
Cov_predict= Mu_x_subtract * diag(Weight) * Mu_x_subtract'+Q;
%Measurement Update
SPs2 = Find_SigmaPoints(Mean_predict, Cov_predict);%sigma points     (n X 2n)
SigmaPoints_Obs = zeros(3, NumOfSP);
for i = 1:NumOfSP
    SigmaPoints_Obs(:, i) = h(SPs2(:, i));
end
Weighted_y = SigmaPoints_Obs* diag(Weight);
Mu_y = sum(Weighted_y, 2); %Z(k ,k-1)  sum(A,2) 求A矩阵的行总和  就是A 的每一行 所有元素求和
Mu_y_subtract = SigmaPoints_Obs - repmat(Mu_y, [1, NumOfSP]);

P_yy = Mu_y_subtract * diag(Weight) * Mu_y_subtract' + R;
P_xy = Mu_x_subtract * diag(Weight) * Mu_y_subtract';

K = P_xy/(P_yy);
Mean_Update = Mean_predict + K*(Observation - Mu_y);%Mu_x: Mean u_x
Cov_Update = Cov_predict-K*P_xy';
end

function [XinPolar] = h(XinCartesian)
x = XinCartesian(1); %第一行那一个数 position x
y = XinCartesian(2);  %第二行那个数  position y
vx= XinCartesian(3);
vy= XinCartesian(4);
XinPolar = [sqrt(x^2+y^2); atan2(y, x);sqrt(vx^2+vy^2)]; % 
end

function [SPs] = Find_SigmaPoints(XinCartesian, Cov)
NumOfDim = size(XinCartesian, 1);%维度是n  % Evaluation points (nx2n)
pnts = [eye(NumOfDim) -eye(NumOfDim)];%n行  2n列 pnts  points
  % Scaling
pnts = sqrt(NumOfDim)*pnts;
pnts= chol(Cov)'*pnts;
SigmaPoints = repmat(XinCartesian, [1,NumOfDim*2]) + pnts;%n*2n  行列数
SPs = SigmaPoints; 
end