T = 1;
F = [1 0 T 0;...
     0 1 0 T;...
     0 0 1 0;...
     0 0 0 1];

%% Target
store_Target_real_state = [];
Target_state = [1980; 2000; 0.2; 0.3];
Target_acceleration = [0; 0; 0; 0];
%% Observer 1
store_Observer_state = [];
store_Target_real_ori = [];
store_Target_real_range = [];
store_Target_real_velocity = [];
Observer_state = [0; 0; 0; 0];
%% Generate Data
% Target
store_Target_real_state = [store_Target_real_state, Target_state];
% Observer 1
store_Observer_state = [store_Observer_state, Observer_state];
relative_pos = [Target_state(1)-Observer_state(1); Target_state(2)-Observer_state(2)];
relative_velocity=[Target_state(3)-Observer_state(3); Target_state(4)-Observer_state(4)];
% Orientation
Target_real_ori = atan2( relative_pos(2), relative_pos(1) );
store_Target_real_ori = [store_Target_real_ori, Target_real_ori];
% Range
Target_real_range = norm(relative_pos);
store_Target_real_range = [store_Target_real_range, Target_real_range];
%Velocity
Target_real_velocity = norm(relative_velocity);
store_Target_real_velocity= [store_Target_real_velocity, Target_real_velocity];
ProcessNoiseCov = 1e-5 * [...
    1/3*T^3, 0,           1/2*T^2,   0;...
    0,       1/3*T^3,     0,         1/2*T^2;...
    1/2*T^2, 0,           T,         0;...
    0,       1/2*T^2,     0,         T];
%go straight
for t = 1:1999
    randPN = mvnrnd([0, 0, 0, 0], ProcessNoiseCov)';
    % Target
    Target_state = F*Target_state + randPN;
    store_Target_real_state = [store_Target_real_state, Target_state];
    % Observer
    Observer_state = F*Observer_state;
    store_Observer_state = [store_Observer_state, Observer_state];
    % Relative Position and Velocity
    relative_pos = [Target_state(1)-Observer_state(1); Target_state(2)-Observer_state(2)];  
    relative_velocity=[Target_state(3)-Observer_state(3); Target_state(4)-Observer_state(4)];
    % Orientation
    Target_real_ori = atan2( relative_pos(2), relative_pos(1) );
    store_Target_real_ori = [store_Target_real_ori, Target_real_ori];  
    % Range
    Target_real_range = norm(relative_pos);
    store_Target_real_range = [store_Target_real_range, Target_real_range];   
     % Velocity
    Target_real_velocity = sqrt((Target_state(3)-Observer_state(3))^2+(Target_state(4)-Observer_state(4))^2);
    store_Target_real_velocity = [store_Target_real_velocity, Target_real_velocity]; 
end
%%
sd_ori = 2e-2; %0.02
sd_range = 1e-3; %0.003
sd_velocity=3e-4; %0.0004

Target_ori = [];
Target_range = [];
Target_velocity=[];
for i = 1:size(store_Target_real_state, 2)
    %relative_pos = [store_Target_state(1, i)-store_Observer_state(1, i); store_Target_state(2, i)-store_Observer_state(2, i)];
    Target_ori = [Target_ori, normrnd(store_Target_real_ori(1, i), sd_ori)];
    Target_range = [Target_range, normrnd(store_Target_real_range(1, i), sd_range)];
    Target_velocity=[Target_velocity,normrnd(store_Target_real_velocity(1,i),sd_velocity)];
end
Observations = [Target_range; Target_ori;Target_velocity];