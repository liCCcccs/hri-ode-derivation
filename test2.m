m1 = 1;
m2 = 2;
K = 10;

A = [0, 0, 0, 1, -1/K;
     0, 0, 1, 0, 1/K;
     0, 0, 0, 0, 1/m1;
     0, 0, 0, 0, -1/m2;
     0, 0, -1/K, 1/K, 0];

B = [0, 0;
     0, 0,;
     1/m1, 0;
     0, 1/m2;
     0, 0];

C = [1, 0, 0, 0, 0;
     0, 0, 0, 0, 1];
D = [0, 0;
     0, 0];


% A = [1.1269   -0.4940    0.1129 
%      1.0000         0         0 
%           0    1.0000         0];
% B = [-0.3832
%       0.5919
%       0.5191];
% C = [1 0 0];
Ob = obsv( A, C);

%% Kalman Filter
Ts = -1;
sys = ss(A,[B B],C,D,Ts,'InputName',{'u1' 'u2' 'w1' 'w2'},'OutputName',{'y1', 'y2'});  % Plant dynamics and additive input noise w

Q = 2.3; 
R = 1; 

[kalmf,L,~,Mx,Z] = kalman(sys,Q,R);
kalmf = kalmf(1,:);

%% Simulate System
sys.InputName = {'u','w'};
sys.OutputName = {'yt1', 'yt1'};
vIn = sumblk({'y1=yt+v'});

kalmf.InputName = {'u','y'};
kalmf.OutputName = 'ye';

SimModel = connect(sys,vIn,kalmf,{'u','w','v'},{'yt','ye'});

t = (0:100)';
u1 = sin(t/5);
u2 = 0*sin(t/5);

rng(10,'twister');
w1 = sqrt(Q)*randn(length(t),1);
w2 = sqrt(Q)*randn(length(t),1);
v = sqrt(R)*randn(length(t),1);

out = lsim(SimModel,[u,w,v]);

yt = out(:,1);   % true response
ye = out(:,2);  % filtered response
y = yt + v;     % measured response

%% Visualise
clf
figure('Renderer', 'painters', 'Position', [300 300 500 500])
subplot(211), plot(t,yt,'b',t,ye,'r--'), 
xlabel('Number of Samples'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered')
subplot(212), plot(t,yt-y,'g',t,yt-ye,'r--'),
xlabel('Number of Samples'), ylabel('Error')
legend('True - measured','True - filtered')





