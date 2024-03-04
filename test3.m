m1 = 1;
m2 = 2;
K = 10;
m1_wrong = 1;
m2_wrong = 4;
K_wrong = 13;

%% System nominal model
A_n = [0, 1;
       0, 0];
B_n = [0, 0;
       -1/m2, 1/m2];
C_n = [1, 0];
D_n = [-1/K, 0];

%% Model in KF
A = [0, 1;
     0, 0];
B = [0, 0;
     -1/m2_wrong, 1/m2_wrong];
C = [1, 0];
D = [-1/K_wrong, 0];

Ob = obsv( A, C);

%% Kalman Filter
Ts = -1;
sys = ss(A_n,[B_n B_n],C_n,[D_n D_n],Ts,'InputName',{'u1' 'u2' 'w1' 'w2'},'OutputName',{'x1', 'x2', 'x3', 'x4', 'x5'});  % Plant dynamics and additive input noise w
sysKF = ss(A,[B B],C,[D D],Ts,'InputName',{'u1' 'u2' 'w1' 'w2'},'OutputName',{'y1', 'y2'});  % Plant dynamics and additive input noise w

Q = [1 0; 0 1]; 
R = [0.001 0; 0 0.001];

[kalmf,L,~,Mx,Z] = kalman(sysKF,Q,R);
%kalmf = kalmf([1,2],:);

%% Simulate System
sys.InputName = {'u1', 'u2', 'w1', 'w2'};
sys.OutputName = {'xt1', 'xt2', 'xt3', 'xt4', 'xt5'};
vIn1 = sumblk('y1 = xt1 + v1');
vIn2 = sumblk('y2 = xt5 + v2');

kalmf.InputName = {'u1', 'u2','y1', 'y2'};
kalmf.OutputName = {'ye1', 'ye2', 'xe1', 'xe2', 'xe3', 'xe4', 'xe5'};

SimModel = connect(sys,vIn1, vIn2, kalmf,{'u1', 'u2', 'w1', 'w2','v1', 'v2'},{'xt1', 'xt2', 'xt3', 'xt4', 'xt5', 'ye1', 'ye2',  'xe1', 'xe2', 'xe3', 'xe4', 'xe5'});

t = (0:100)';
u1 = 10*sin(t/10);
u2 = 0*sin(t);

rng(10,'twister');
w1 = zeros(length(t),1); %sqrt(Q(1,1))*randn(length(t),1);
w2 = zeros(length(t),1); %sqrt(Q(2,2))*randn(length(t),1);
v1 = zeros(length(t),1); %sqrt(R(1,1))*randn(length(t),1);
v2 = zeros(length(t),1); %sqrt(R(2,2))*randn(length(t),1);

out = lsim(SimModel,[u1, u2, w1, w2, v1, v2]);

xt1 = out(:,1);   % true response
xt2 = out(:,2);   % true response
xt3 = out(:,3);   % true response
xt4 = out(:,4);   % true response
xt5 = out(:,5);   % true response

ye1 = out(:,6);  % filtered response
ye2 = out(:,7);  % filtered response
xe1 = out(:,8);  % filtered response
xe2 = out(:,9);  % filtered response
xe3 = out(:,10);  % filtered response
xe4 = out(:,11);  % filtered response
xe5 = out(:,12);  % filtered response

y1 = xt1 + v1;     % measured response
y2 = xt2 + v2;     % measured response

%% Visualise
clf
figure('Renderer', 'painters', 'Position', [300 300 800 500])
plot(t,xt3,'b',t,xe3,'r--'), 
xlabel('Number of Samples'), ylabel('Output')
title('Kalman Filter Response')
legend('True','Filtered')
% subplot(212), plot(t,yt-y,'g',t,yt-ye,'r--'),
% xlabel('Number of Samples'), ylabel('Error')
% legend('True - measured','True - filtered')

%%
A = [0, 1;
    0, 0];
C = [1, 0];
Ob = obsv(A,C)



