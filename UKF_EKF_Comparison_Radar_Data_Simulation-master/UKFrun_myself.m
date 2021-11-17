function [ Mean_Update, Cov_Update ] = UKFrun_myself( Mean, Cov, Observation, Q, R, T )
F = [1, 0, T, 0;
    0, 1, 0, T;
    0, 0, 1, 0;
    0, 0, 0, 1];
% Parameters for UKF
n = size(Mean, 1);
Alpha =0.25;
K = 2;
Beta = 2;%4
Lambda = Alpha^2*(n+K)-n;
C = sqrt(Lambda+n);
SPs = Find_SigmaPoints(Mean, Cov, C );%sigma points
NumOfSP = size(SPs, 2);%number of sigma points
%      Mean_predict=zeros(n,1);
%weights
W_Mean = [Lambda/(2*n+Lambda), repmat(1/(2*(2*n+Lambda)), [1, 2*n])];% 1*(2n+1) 
W_Cov = [Lambda/(2*n+Lambda)+(1-Alpha^2+Beta), repmat(1/(2*(n+Lambda)), [1, 2*n])];
%       Cov_predict=zeros(n,n);
Mean_predict_star=F*SPs;

Weighted_x = Mean_predict_star * diag(W_Mean);
Mean_predict = sum(Weighted_x, 2);
Mu_x_subtract = Mean_predict_star - repmat(Mean_predict, [1, NumOfSP]);
Cov_predict= Mu_x_subtract * diag(W_Cov) * Mu_x_subtract'+Q;

SigmaPoints_Obs = zeros(3, NumOfSP);
for i = 1:NumOfSP
    SigmaPoints_Obs(:, i) = h(Mean_predict_star(:, i));
end
Weighted_y = SigmaPoints_Obs* diag(W_Mean);
Mu_y = sum(Weighted_y, 2); %Z(k ,k-1)  sum(A,2) ��A��������ܺ�  ����A ��ÿһ�� ����Ԫ�����
Mu_y_subtract = SigmaPoints_Obs - repmat(Mu_y, [1, NumOfSP]);

P_yy = Mu_y_subtract * diag(W_Cov) * Mu_y_subtract' + R;
P_xy = Mu_x_subtract * diag(W_Cov) * Mu_y_subtract';

K = P_xy*inv(P_yy);
Mean_Update = Mean_predict + K*(Observation - Mu_y);%Mu_x: Mean u_x
Cov_Update = Cov_predict-P_xy*inv(P_yy)*P_xy';% !!!!!�ر���� ���鱾�Ϲ�ʽ��ȫ��һ�� �����Ϲ�ʽ�ͽ�������� ����ԭ�򣡣�
%����һ���������ҵ��˸ñ����ʽ��update��ʽ�����շ���
end
function [XinPolar] = h(XinCartesian)
x = XinCartesian(1); %��һ����һ���� position x
y = XinCartesian(2);  %�ڶ����Ǹ���  position y
vx= XinCartesian(3);
vy= XinCartesian(4);
XinPolar = [sqrt(x^2+y^2); atan2(y, x);sqrt(vx^2+vy^2)]; % 
end

function [SPs] = Find_SigmaPoints(XinCartesian, Cov, C)
%  Using eigenvalue decomposition ����ֵ�ֽ⣬����ֵ�ֽ⣨SVD��
NumOfDim = size(XinCartesian, 1);%ά����n
NumOfPoints = NumOfDim*2;%2n
[V, D] = eig(Cov);
pnts = zeros(NumOfDim, NumOfPoints);%n��  2n�� pnts  points
for i = 1:NumOfDim
    pnts(i, i) = sqrt(abs(D(i, i)));
    pnts(i, i+NumOfDim) = -sqrt(abs(D(i, i)));
end
%scale samples such that they have the same second order states as a standard normal
Covariance_pnts = C*pnts;
Rotated_pnts = V*Covariance_pnts;
SigmaPoints = repmat(XinCartesian, [1, NumOfPoints]) + Rotated_pnts;%n*2n  ������
SPs = [XinCartesian, SigmaPoints]; %���Բ�û��������ʼֵ�� ������Ŀȷʵ��2n+1��
end