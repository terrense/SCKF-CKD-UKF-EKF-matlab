function [ Mean_Update, Cov_Update,S_Update] = SCKFrun( Mean, ~, Observation, Q, R, T,S_input )
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
S=S_input;
SPs1 = Find_SigmaPoints(Mean, S);%sigma points     (n X 2n)
NumOfSP = size(SPs1, 2);%number of sigma points
Mean_predict_star=F*SPs1;
Weighted_x = Mean_predict_star * diag(Weight);
Mean_predict = sum(Weighted_x, 2);
Mu_x_subtract = Mean_predict_star - repmat(Mean_predict, [1, NumOfSP]);
ka_fang_star=(1/(sqrt(NumOfSP)))*Mu_x_subtract;
SQ=Find_SQ_SR(Q);
[~,l]=qr([ka_fang_star,SQ]);

%Measurement Update
SPs2 = Find_SigmaPoints(Mean_predict, l);%sigma points     (n X 2n)
SigmaPoints_Obs = zeros(2, NumOfSP);
for i = 1:NumOfSP
    SigmaPoints_Obs(:, i) = h(SPs2(:, i));
end
Weighted_y = SigmaPoints_Obs* diag(Weight);
Mu_y = sum(Weighted_y, 2); %Z(k ,k-1)  sum(A,2) 求A矩阵的行总和  就是A 的每一行 所有元素求和
ga_ma =(1/(sqrt(NumOfSP)))*SigmaPoints_Obs - repmat(Mu_y, [1, NumOfSP]);
ka_fang=(1/(sqrt(NumOfSP)))*SPs2 - repmat(Mean_predict, [1, NumOfSP]);
SR=Find_SQ_SR(R);
P_xy = ka_fang * ga_ma';
[~,S_yy]=qr([ga_ma,SR]);
% P_yy = S_yy'*S_yy;
K = (P_xy/S_yy')/S_yy;
% K = (mldivide(P_yy',P_xy'))';
% K = mrdivide(P_xy,P_yy);
Mean_Update = Mean_predict + K*(Observation - Mu_y);%Mu_x: Mean u_x
[~,S_Update]=qr([ka_fang-K*ga_ma, K*SR]);
Cov_Update = S_Update'*S_Update;
end

function [XinPolar] = h(XinCartesian)
x = XinCartesian(1); %第一行那一个数 position x
y = XinCartesian(2);  %第二行那个数  position y
XinPolar = [sqrt(x^2+y^2); atan2(y, x)]; % 距离 还有 方向角
end

function [SPs] = Find_SigmaPoints(XinCartesian, S)
NumOfDim = size(XinCartesian, 1);%维度是n
NumOfPoints = NumOfDim*2;%2n
pnts = zeros(NumOfDim, NumOfPoints);%n行  2n列 pnts  points
for i = 1:NumOfDim
    pnts(i, i) = abs(S(i, i));
    pnts(i, i+NumOfDim) = -abs(S(i, i));
end
Covariance_pnts = pnts;
Rotated_pnts =Covariance_pnts;
SigmaPoints = repmat(XinCartesian, [1, NumOfPoints]) + Rotated_pnts;%n*2n  行列数
SPs = SigmaPoints; 
end
function [SQ_R] = Find_SQ_SR(Cov)
NumOfDim = size(Cov, 1);%维度是n
pnts = zeros(NumOfDim, NumOfDim);%n行  2n列 pnts  points
for i = 1:NumOfDim
  for j = 1:NumOfDim
      pnts(i, j) = sqrt(abs(Cov(i, j)));
  end
end
SQ_R =pnts; 
end