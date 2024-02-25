syms l1 l2 l3 l4
syms lc1 la1 lb1 lc2 lc3 lc4 la4 la2
syms q1 h_q2 dq1 h_dq2 ddq1 h_ddq2 
syms r_d2 r_dd2 r_ddd2
syms r_d3 r_dd3 r_ddd3
syms r_q4 r_dq4 r_ddq4
syms r_q5 r_dq5 r_ddq5
syms K_AFz K_AFx K_AMy
syms K_BFz K_BFx K_BMy
syms K_CFz K_CFx K_CMy

syms K_AFzN3 K_AFxN3 K_AMyN3
syms K_BFzN3 K_BFxN3 K_BMyN3

syms D_AFz D_AFx D_AMy
syms D_BFz D_BFx D_BMy
syms D_CFz D_CFx D_CMy
syms m1 m2 m3 m4 g
syms I_G1z I_G2z I_G3z I_G4z
syms tau1 tau2 tau3 tau4

q1 = pi/2; dq1 = 0; ddq1 = 0;
tau1 = 0; tau3 = 0;
 

% intermidiate variables
syms f_intBi  % D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))))
syms f_intBj  % - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))))
syms f_intCi  % D_CFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_CFx*((l2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (l2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + l4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + l4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))))
syms f_intCj  % - D_CFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_CFz*((l2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (l2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + l4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + l4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))))

syms m4gi  % m4*g*cos(q1+r_q4+r_q5)
syms m4gj  % -m4*g*sin(q1+r_q4+r_q5)
syms tau_intBz  % K_BMy*(r_q4 + r_q5 - h_q2) + D_BMy*(r_dq4 + r_dq5 - h_dq2)
syms tau_intCz  % D_CMy*(r_dq4 - h_dq2 + r_dq5) + K_CMy*(r_q4 - h_q2 + r_q5)

syms f_intAi  % K_AFx * r_d3 + D_AFx * r_dd3
syms f_intAj  % + K_AFz * r_d2 + D_AFz * r_dd2
syms m3gi  % m3*g*cos(q1+r_q4)
syms m3gj  % -m3*g*sin(q1+r_q4)
syms tau_intAz  % K_AMy*r_q4 + D_AMy*r_dq4
syms m2gi  % m2*g*cos(q1+h_q2)
syms m2gj  % -m2*g*sin(q1+h_q2)
syms m1gi  % m1*g*cos(q1)
syms m1gj  % -m1*g*sin(q1)

% test
syms a1 a2
syms b1 b2 b3 b4 b5
syms c1 c2 c3 c4 c5
syms d1 d2 d3 d5
syms e1 e2 e3 e4 e5

eqns = [
    tau2 == a1 + a2*h_ddq2

    0 == b1 + r_ddd3*(b3) + r_ddd2*(b2) + r_ddq4*(b4) + r_ddq5*(b5)  % no active input

    tau4 == c1 + r_ddq5*(c5) + r_ddq4*(c4) + r_ddd3*c3 + r_ddd2*c2

    0 == d1 + (d2)*r_ddd2 + (d3)*r_ddd3 + (d5)*r_ddq5

    0 == e1 + r_ddd3*(e3) + r_ddd2*(e2) + r_ddq4*(e4) + r_ddq5*e5

    ];

% ddthetaSolved = solve(eqns, [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5]);


% simplify(ddthetaSolved.h_ddq2)
% simplify(ddthetaSolved.r_ddd2)
% simplify(ddthetaSolved.r_ddd3)
% simplify(ddthetaSolved.r_ddq4)
% simplify(ddthetaSolved.r_ddq5)


vSym.m1 = 7.275;
vSym.m2 = 3.75;
vSym.m3 = 2;
vSym.m4 = 2;
vSym.g = 10;
vSym.I_G1z = 0.121; vSym.I_G2z = 0.055; vSym.I_G3z = 0.02; vSym.I_G4z = 0.02; 
vSym.l1 = 0.4; vSym.l2 = 0.4; vSym.l3 = 0.2; vSym.l4 = 0.4;
vSym.lc1 = 0.173; vSym.la1 = 0.2; vSym.lb1 = 0.2; vSym.lc2 = 0.173; 
vSym.lc3 = 0.1; vSym.lc4 = 0.2; vSym.la4 = 0.2; vSym.la2 = 0.2; 
vSym.K_AFz = 6000; vSym.K_AFx = 3000; vSym.K_AMy = 50;
vSym.K_BFz = 6000; vSym.K_BFx = 3000; vSym.K_BMy = 500;

vSym.K_AFzN3 = 60000*5; vSym.K_AFxN3 = 30000*5; vSym.K_AMyN3 = 100;
vSym.K_BFzN3 = 60000*5; vSym.K_BFxN3 = 30000*5; vSym.K_BMyN3 = 100;

vSym.K_CFz = 24000; vSym.K_CFx = 12000; vSym.K_CMy = 200;
vSym.D_AFz = 100; vSym.D_AFx = 100; vSym.D_AMy = 10;
vSym.D_BFz = 100; vSym.D_BFx = 100; vSym.D_BMy = 10;
vSym.D_CFz = 100; vSym.D_CFx = 100; vSym.D_CMy = 10;


% Wrong model parameter
vSym_est.m1 = 7.275;
vSym_est.m2 = 3.75;
vSym_est.m3 = 2;
vSym_est.m4 = 2;
vSym_est.g = 10;
vSym_est.I_G1z = 0.121; vSym_est.I_G2z = 0.055; vSym_est.I_G3z = 0.02; vSym_est.I_G4z = 0.02; 
vSym_est.l1 = 0.4; vSym_est.l2 = 0.4; vSym_est.l3 = 0.2; vSym_est.l4 = 0.4;
vSym_est.lc1 = 0.173; vSym_est.la1 = 0.2; vSym_est.lb1 = 0.2; vSym_est.lc2 = 0.173; 
vSym_est.lc3 = 0.1; vSym_est.lc4 = 0.2; vSym_est.la4 = 0.2; vSym_est.la2 = 0.2; 
vSym_est.K_AFz = 3000; vSym_est.K_AFx = 2000; vSym_est.K_AMy = 50;
vSym_est.K_BFz = 3000; vSym_est.K_BFx = 2000; vSym_est.K_BMy = 500;

vSym_est.K_AFzN3 = 30000*5; vSym_est.K_AFxN3 = 20000*5; vSym_est.K_AMyN3 = 100;
vSym_est.K_BFzN3 = 30000*5; vSym_est.K_BFxN3 = 20000*5; vSym_est.K_BMyN3 = 100;

vSym_est.K_CFz = 24000; vSym_est.K_CFx = 12000; vSym_est.K_CMy = 200;
vSym_est.D_AFz = 100; vSym_est.D_AFx = 100; vSym_est.D_AMy = 10;
vSym_est.D_BFz = 100; vSym_est.D_BFx = 100; vSym_est.D_BMy = 10;
vSym_est.D_CFz = 100; vSym_est.D_CFx = 100; vSym_est.D_CMy = 10;

%% test
% t = 0;
% y0 = [q1, dq1, -pi/2, 0, 0, 0, 0, 0, 0, 0, -pi/2, 0];  %initial condition
% tau = [0, 0, 0, 0];
% dydt = runrobot(t, y0, tau, vSym)

%% visualise 
% x = 0:0.001:0.1;
% ylin = vSym.K_AFz * x;
% ynonlin = vSym.K_AFz * x + vSym.K_AFzN3*5 * x.^3; 
% figure('Renderer', 'painters', 'Position', [300 300 500 500])
% hold on
% plot(x*1000, ylin)
% plot(x*1000, ynonlin)
% xlabel("Displacement (mm)"); ylabel("Force (N)");

%% Calculate equilibrium starting posiiton
syms x_eq; % Define x_eq as symbolic variable
equation = vSym.K_AFz * x_eq + vSym.K_BFx * x_eq == (vSym.m3 + vSym.m4) * vSym.g;
x_eq_val = double(solve(equation, x_eq));
%% solve the dynamic equations
tspan = [0 :0.04: 10];
%tspan = 0:0.05:1;
y0 = [q1, dq1, -pi/2, 0, 0, 0, x_eq_val, 0, 0, 0, -pi/2, 0];  %initial condition
tau = [0, 0, 0, 0];
y = zeros(length(tspan), length(y0));
y(1, :) = y0;
for i = 2:length(tspan)
    ts = [0, 0.04];
    y0 = y(i-1, :);
    [t_tmp,y_tmp] = ode23(@(t_tmp,y_tmp) runrobot(t_tmp,y_tmp, tau, vSym), ts, y0);
    y(i, :) = y_tmp(end, :);
end

% tspan = [0 :0.04: 2];
% % tspan = 0:0.05:1;
% y0 = [q1, dq1, -pi/2, 0, 0, 0, x_eq_val, 0, 0, 0, -pi/2, 0];  %initial condition
% tau = [0, 0, 0, 0];
% [t,y] = ode23(@(t,y) runrobot(t,y, tau, vSym), tspan, y0);

%% observer
tspan = [0 :0.04: 10];
%tspan = 0:0.05:1;
y0 = [q1, dq1, -pi/2, 0, 0, 0, x_eq_val, 0, 0, 0, -pi/2, 0];  %initial condition
tau = [0, 0, 0, 0];
y_ob = zeros(length(tspan), length(y0));
y_ob(1, :) = y0;
%tspan = 0:0.05:1;
for i = 2:length(tspan)
    ts = [0, 0.04];
    y0 = y_ob(i-1, :);
    y_true = y(i, :);
    [t_tmp,y_tmp] = ode23(@(t_tmp,y_tmp) runobserver(t_tmp,y_tmp, tau, vSym_est, y_true), ts, y0);
    y_ob(i, :) = y_tmp(end, :);
end

% y0 = [q1, dq1, -pi/2, 0, 0, 0, x_eq_val, 0, 0, 0, -pi/2, 0];  %initial condition
% tau = [0, 0, 0, 0];
% [t_ob,y_ob] = ode23(@(t,y) runrobot(t,y, tau, vSym_est), tspan, y0);

%% Plot trajectory
% figure('Renderer', 'painters', 'Position', [300 300 800 800])
% plot(t, y(:, [3, 11]))
% legend(["h th", "r th"])
% % 
% q1 = y(:, 1); h_q2 = y(:, 3); r_d2 = y(:, 5); r_d3 = y(:, 7); r_q4 = y(:, 9); r_q5 = y(:, 11);
% dq1 = y(:, 2); h_dq2 = y(:, 4); r_dd2 = y(:, 6); r_dd3 = y(:, 8); r_dq4 = y(:, 10); r_dq5 = y(:, 12);
% 
% l1 = vSym.l1;  l2 = vSym.l2;  l3 = vSym.l3;  l4 = vSym.l4;
% lc1 = vSym.lc1;  la1 = vSym.la1;  lb1 = vSym.lb1;  lc2 = vSym.lc2;
% lc3 = vSym.lc3;  lc4 = vSym.lc4;  la4 = vSym.la4;  la2 = vSym.la2;
% K_AFz = vSym.K_AFz;  K_AFx = vSym.K_AFx;  K_AMy = vSym.K_AMy;
% K_BFz = vSym.K_BFz;  K_BFx = vSym.K_BFx;  K_BMy = vSym.K_BMy;
% D_AFz = vSym.D_AFz;  D_AFx = vSym.D_AFx;  D_AMy = vSym.D_AMy;
% D_BFz = vSym.D_BFz;  D_BFx = vSym.D_BFx;  D_BMy = vSym.D_BMy;
% m1 = vSym.m1;  m2 = vSym.m2;  m3 = vSym.m3;  m4 = vSym.m4;  g = vSym.g;
% I_G1z = vSym.I_G1z;  I_G2z = vSym.I_G2z;  I_G3z = vSym.I_G3z;  I_G4z = vSym.I_G4z;
% 
% f_intBi = D_BFx.*(r_dd2.*cos(h_q2) + r_dd3.*cos(h_q2) - h_dq2.*l1.*sin(h_q2) + r_d3.*r_dq4.*cos(h_q2) + r_d3.*r_dq5.*cos(h_q2) + la1.*r_dq4.*sin(h_q2) + la1.*r_dq5.*sin(h_q2) + r_d2.*r_dq4.*sin(h_q2) + r_d2.*r_dq5.*sin(h_q2) + l3.*r_dq5.*sin(h_q2 - r_q4)) - K_BFx.*((la2.*(cos(h_q2).*cos(q1 - pi/2) - sin(h_q2).*sin(q1 - pi/2)) + l1.*cos(q1 - pi/2)).*(cos(h_q2).*cos(q1 - pi/2) - sin(h_q2).*sin(q1 - pi/2)) + (la2.*(cos(h_q2).*sin(q1 - pi/2) + sin(h_q2).*cos(q1 - pi/2)) + l1.*sin(q1 - pi/2)).*(cos(h_q2).*sin(q1 - pi/2) + sin(h_q2).*cos(q1 - pi/2)) + (cos(h_q2).*sin(q1 - pi/2) + sin(h_q2).*cos(q1 - pi/2)).*(l3.*(cos(q1).*cos(r_q4) - sin(q1).*sin(r_q4)) + cos(q1).*(la1 + r_d2) + r_d3.*sin(q1) + la4.*(cos(r_q5).*(cos(q1).*cos(r_q4) - sin(q1).*sin(r_q4)) - sin(r_q5).*(cos(q1).*sin(r_q4) + cos(r_q4).*sin(q1)))) - (cos(h_q2).*cos(q1 - pi/2) - sin(h_q2).*sin(q1 - pi/2)).*(l3.*(cos(q1).*sin(r_q4) + cos(r_q4).*sin(q1)) + sin(q1).*(la1 + r_d2) - r_d3.*cos(q1) + la4.*(cos(r_q5).*(cos(q1).*sin(r_q4) + cos(r_q4).*sin(q1)) + sin(r_q5).*(cos(q1).*cos(r_q4) - sin(q1).*sin(r_q4)))));
% f_intBj = - D_BFz.*(r_dd2.*sin(h_q2) + r_dd3.*sin(h_q2) - la1.*r_dq4.*cos(h_q2) - la1.*r_dq5.*cos(h_q2) - r_d2.*r_dq4.*cos(h_q2) - r_d2.*r_dq5.*cos(h_q2) + r_d3.*r_dq4.*sin(h_q2) + r_d3.*r_dq5.*sin(h_q2) - l3.*r_dq5.*cos(h_q2 - r_q4) + h_dq2.*l1.*cos(h_q2)) - K_BFz.*((la2.*(cos(h_q2).*sin(q1 - pi/2) + sin(h_q2).*cos(q1 - pi/2)) + l1.*sin(q1 - pi/2)).*(cos(h_q2).*cos(q1 - pi/2) - sin(h_q2).*sin(q1 - pi/2)) - (la2.*(cos(h_q2).*cos(q1 - pi/2) - sin(h_q2).*sin(q1 - pi/2)) + l1.*cos(q1 - pi/2)).*(cos(h_q2).*sin(q1 - pi/2) + sin(h_q2).*cos(q1 - pi/2)) + (cos(h_q2).*sin(q1 - pi/2) + sin(h_q2).*cos(q1 - pi/2)).*(l3.*(cos(q1).*sin(r_q4) + cos(r_q4).*sin(q1)) + sin(q1).*(la1 + r_d2) - r_d3.*cos(q1) + la4.*(cos(r_q5).*(cos(q1).*sin(r_q4) + cos(r_q4).*sin(q1)) + sin(r_q5).*(cos(q1).*cos(r_q4) - sin(q1).*sin(r_q4)))) + (cos(h_q2).*cos(q1 - pi/2) - sin(h_q2).*sin(q1 - pi/2)).*(l3.*(cos(q1).*cos(r_q4) - sin(q1).*sin(r_q4)) + cos(q1).*(la1 + r_d2) + r_d3.*sin(q1) + la4.*(cos(r_q5).*(cos(q1).*cos(r_q4) - sin(q1).*sin(r_q4)) - sin(r_q5).*(cos(q1).*sin(r_q4) + cos(r_q4).*sin(q1)))));
% tau_intBz  = K_BMy.*(r_q4 + r_q5 - h_q2) + D_BMy.*(r_dq4 + r_dq5 - h_dq2);
% f_intAi  = - K_AFx .* r_d2 - D_AFx .* r_dd2;
% f_intAj = - K_AFz .* r_d3 - D_AFz .* r_dd3;
% tau_intAz  = K_AMy.*r_q4 + D_AMy.*r_dq4;
% 
% tau_h_eqv = f_intBj * vSym.la2 + f_intAj * vSym.la4;
% ang_diff = rad2deg(r_q5 - h_q2);
% 
% plot(f_intAj)
% hold on
% plot(f_intBj)
% plot(tau_h_eqv)
% plot(ang_diff)
% legend("f_intBi", "f_intBj", "tau_h_eqv", "ang_diff")
% 
% figure('Renderer', 'painters', 'Position', [300 300 800 800])
% scatter(ang_diff, tau_h_eqv)

%% save data
% hri_data.ang_diff = rad2deg(r_q5 - h_q2);
% hri_data.vel_diff = rad2deg(r_dq5 - h_dq2);
% hri_data.r_q5 = r_q5;
% hri_data.h_q2 = h_q2;
% hri_data.tau_h_eqv = tau_h_eqv;
% save('hri_data.mat', 'hri_data')

%% Visualise
figure('Renderer', 'painters', 'Position', [300 300 800 800])

dt = 0.04;
q1 = y(:, 1); h_q2 = y(:, 3); r_d2 = y(:, 5); r_d3 = y(:, 7); r_q4 = y(:, 9); r_q5 = y(:, 11);
dq1 = y(:, 2); h_dq2 = y(:, 4); r_dd2 = y(:, 6); r_dd3 = y(:, 8); r_dq4 = y(:, 10); r_dq5 = y(:, 12);

q1_est = y_ob(:, 1); h_q2_est = y_ob(:, 3); r_d2_est = y_ob(:, 5); r_d3_est = y_ob(:, 7); r_q4_est = y_ob(:, 9); r_q5_est = y_ob(:, 11);
dq1_est = y_ob(:, 2); h_dq2_est = y_ob(:, 4); r_dd2_est = y_ob(:, 6); r_dd3_est = y_ob(:, 8); r_dq4_est = y_ob(:, 10); r_dq5_est = y_ob(:, 12);

symbols = {m1, m2, m3, m4, g, I_G1z, I_G2z, I_G3z, I_G4z, ...
           l1, l2, l3, l4,...
           lc1, la1, lb1, lc2, lc3, lc4, la4, la2,...
           K_AFz, K_AFx, K_AMy,...
           K_BFz, K_BFx, K_BMy, ...
           D_AFz, D_AFx, D_AMy,...
           D_BFz, D_BFx, D_BMy, ...
           q1, dq1, ...
           h_q2, h_dq2,...
           r_d2, r_dd2,...
           r_d3, r_dd3,...
           r_q4, r_dq4,...
           r_q5, r_dq5,...
           tau1, tau2, tau3, tau4};
    symbol_vals = {vSym.m1, vSym.m2, vSym.m3, vSym.m4, vSym.g,...
                   vSym.I_G1z, vSym.I_G2z, vSym.I_G3z, vSym.I_G4z, ...
                   vSym.l1, vSym.l2, vSym.l3, vSym.l4,...
                   vSym.lc1, vSym.la1, vSym.lb1, vSym.lc2, ...
                   vSym.lc3, vSym.lc4, vSym.la4, vSym.la2, ...
                   vSym.K_AFz, vSym.K_AFx, vSym.K_AMy,...
                   vSym.K_BFz, vSym.K_BFx, vSym.K_BMy, ...
                   vSym.D_AFz, vSym.D_AFx, vSym.D_AMy,...
                   vSym.D_BFz, vSym.D_BFx, vSym.D_BMy};
    variable_vals = {q1, dq1, ...
                     h_q2, h_dq2,...
                     r_d2, r_dd2,...
                     r_d3, r_dd3,...
                     r_q4, r_dq4,...
                     r_q5, r_dq5,...
                     tau(1), tau(1), tau(3), tau(4)};

f_intB = zeros(2, length(y));
%f_intB(1, i) = subs(D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4))))), symbols, [symbol_vals, variable_vals]);
%f_intB(1, i) = subs(- D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1))))), symbols, [symbol_vals, variable_vals]);
    

%figure('Renderer', 'painters', 'Position', [300 300 600 500])
%v = VideoWriter('Against10Nm_Neg.avi', 'Motion JPEG AVI');
%v.FrameRate = 25;
%open(v);
for i = 1:length(y)
    [j1, h_j2, h_j3, r_j2, r_j3, r_j4, r_j5, r_j6, h_cuff_A, h_cuff_B, h_cuff_C, r_cuff_B, T] = calc_joint_position(vSym, q1(i), h_q2(i), r_d2(i), r_d3(i), r_q4(i), r_q5(i));

    [j1_est, h_j2_est, h_j3_est, r_j2_est, r_j3_est, r_j4_est, r_j5_est, r_j6_est, h_cuff_A_est, h_cuff_B_est, h_cuff_C_est, r_cuff_B_est, T_est] = calc_joint_position(vSym_est, q1_est(i), h_q2_est(i), r_d2_est(i), r_d3_est(i), r_q4_est(i), r_q5_est(i));


    clf
    plot_2d_links({j1, h_j2, h_j3}, '-');
    plot_2d_links({r_j3, r_j4, r_j5, r_j6}, '-');
    plot_2d_links({h_cuff_A, r_j3}, ':');
    plot_2d_links({h_cuff_B, r_cuff_B}, ':');
    plot_2d_links({h_cuff_C, r_j6}, ':');
    colors = ['b', 'b', 'b', 'r', 'r', 'r'];
    plot_point({h_cuff_A, h_cuff_B, h_cuff_C, r_j3, r_cuff_B, r_j6}, colors)

    plot_2d_links({j1_est, h_j2_est, h_j3_est}, '--');
    plot_2d_links({r_j3_est, r_j4_est, r_j5_est, r_j6_est}, '--');
    plot_2d_links({h_cuff_A_est, r_j3_est}, ':');
    plot_2d_links({h_cuff_B_est, r_cuff_B_est}, ':');
    plot_2d_links({h_cuff_C_est, r_j6_est}, ':');


    view(0, 90)
    legend("Human", "Robot")

    axis([-0.2,0.8,-0.6,0.4, -1,1])
    %frame = getframe(gcf); 
    %writeVideo(v, frame);
    pause(dt)
end
close(v);


%figure('Renderer', 'painters', 'Position', [300 300 800 700])
%plot(f_intB)

%% Calculate differetial equations
function dydt = runrobot(t, y, tau, vSym)

    vVar.q1 = y(1); vVar.dq1 = y(2); 
    vVar.h_q2 = y(3); vVar.h_dq2 = y(4);
    vVar.r_d2 = y(5); vVar.r_dd2 = y(6);
    vVar.r_d3 = y(7); vVar.r_dd3 = y(8);
    vVar.r_q4 = y(9); vVar.r_dq4 = y(10);
    vVar.r_q5 = y(11); vVar.r_dq5 = y(12); 

    vVar.tau1 = tau(1); vVar.tau3 = tau(3);
    if t < 4
        vVar.tau2 = -20;
        vVar.tau4 = 23;
    else
        vVar.tau2 = 0;
        vVar.tau4 = 0;
    end

    %ddq1 = calc_equation_of_motion(ddthetaSolved.ddq1, symbols, vSym, vVar);
    ddq1 = 0;
    [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5] = calc_equation_of_motion(vSym, vVar);

    dydt = zeros(12,1);
    dydt(1) = vVar.dq1; dydt(2) = ddq1;
    dydt(3) = vVar.h_dq2; dydt(4) = h_ddq2;
    dydt(5) = vVar.r_dd2; dydt(6) = r_ddd2;
    dydt(7) = vVar.r_dd3; dydt(8) = r_ddd3;
    dydt(9) = vVar.r_dq4; dydt(10) = r_ddq4;
    dydt(11) = vVar.r_dq5; dydt(12) = r_ddq5;
end

function dydt = runobserver(t, y, tau, vSym, y_true)

    vVar.q1 = y(1); vVar.dq1 = y(2); 
    vVar.h_q2 = y(3); vVar.h_dq2 = y(4);
    vVar.r_d2 = y(5); vVar.r_dd2 = y(6);
    vVar.r_d3 = y(7); vVar.r_dd3 = y(8);
    vVar.r_q4 = y(9); vVar.r_dq4 = y(10);
    vVar.r_q5 = y(11); vVar.r_dq5 = y(12); 

    vVar.tau1 = tau(1); vVar.tau3 = tau(3);
    if t < 4
        %vVar.tau2 = -20;
        vVar.tau2 = -19;
        vVar.tau4 = 23;
    else
        vVar.tau2 = 0;
        vVar.tau4 = 0;
    end

    % % Luenburger observer with 1 measurement (theta_r)
    % L = [1,1,1,1,1]*10;
    % ddq1 = 0;
    % r_q5_true = y_true(11);
    % [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5] = calc_equation_of_motion(vSym, vVar);
    % h_ddq2 = h_ddq2 + L(1) * (r_q5_true - vVar.r_q5);
    % r_ddd2 = r_ddd2 + L(2) * (r_q5_true - vVar.r_q5);
    % r_ddd3 = r_ddd3 + L(3) * (r_q5_true - vVar.r_q5);
    % r_ddq4 = r_ddq4 + L(4) * (r_q5_true - vVar.r_q5);
    % r_ddq5 = r_ddq5 + L(5) * (r_q5_true - vVar.r_q5);

    % Luenburger observer with 7 measurement (theta_r, F_int)
    L = ones(5,7)*0.08;
    L(:, 1) = 10;  % for the robot joint angle
    %L(5, :) = 5;  % for the robot joint angle

    [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5] = calc_equation_of_motion(vSym, vVar);  % estimation
    q1_true = y_true(:, 1); h_q2_true = y_true(:, 3); r_d2_true = y_true(:, 5); r_d3_true = y_true(:, 7); r_q4_true = y_true(:, 9); r_q5_true = y_true(:, 11);
    dq1_true = y_true(:, 2); h_dq2_true = y_true(:, 4); r_dd2_true = y_true(:, 6); r_dd3_true = y_true(:, 8); r_dq4_true = y_true(:, 10); r_dq5_true = y_true(:, 12);

    [f_intAi_true, f_intAj_true, tau_intAz_true, f_intBi_true, f_intBj_true, tau_intBz_true] = calcualte_interaction_force_torque(vSym, q1_true, h_q2_true, r_d2_true, r_d3_true, r_q4_true, r_q5_true, dq1_true, h_dq2_true, r_dd2_true, r_dd3_true, r_dq4_true, r_dq5_true);
    [f_intAi_est, f_intAj_est, tau_intAz_est, f_intBi_est, f_intBj_est, tau_intBz_est] = calcualte_interaction_force_torque(vSym, vVar.q1, vVar.h_q2, vVar.r_d2, vVar.r_d3, vVar.r_q4, vVar.r_q5, vVar.dq1, vVar.h_dq2, vVar.r_dd2, vVar.r_dd3, vVar.r_dq4, vVar.r_dq5);

    error_term = [r_q5_true - vVar.r_q5;
                  f_intAi_true - f_intAi_est;
                  f_intAj_true - f_intAj_est;
                  tau_intAz_true - tau_intAz_est;
                  f_intBi_true - f_intBi_est;
                  f_intBj_true - f_intBj_est;
                  tau_intBz_true - tau_intBz_est];
    lb_term = L*error_term;
    h_ddq2 = h_ddq2 + lb_term(1);
    r_ddd2 = r_ddd2 + lb_term(2);
    r_ddd3 = r_ddd3 + lb_term(3);
    r_ddq4 = r_ddq4 + lb_term(4);
    r_ddq5 = r_ddq5 + lb_term(5);
    ddq1 = 0;

    dydt = zeros(12,1);
    dydt(1) = vVar.dq1; dydt(2) = ddq1;
    dydt(3) = vVar.h_dq2; dydt(4) = h_ddq2;
    dydt(5) = vVar.r_dd2; dydt(6) = r_ddd2;
    dydt(7) = vVar.r_dd3; dydt(8) = r_ddd3;
    dydt(9) = vVar.r_dq4; dydt(10) = r_ddq4;
    dydt(11) = vVar.r_dq5; dydt(12) = r_ddq5;

    % Can we make use of the measured r_q5?
    % dydt(11) = r_q5_true;  % this causes system unstable
end

function [f_intAi, f_intAj, tau_intAz, f_intBi, f_intBj, tau_intBz] = calcualte_interaction_force_torque(vSym, q1, h_q2, r_d2, r_d3, r_q4, r_q5, dq1, h_dq2, r_dd2, r_dd3, r_dq4, r_dq5)
    % Assign values to symbols
    l1 = vSym.l1;  l2 = vSym.l2;  l3 = vSym.l3;  l4 = vSym.l4;
    lc1 = vSym.lc1;  la1 = vSym.la1;  lb1 = vSym.lb1;  lc2 = vSym.lc2;
    lc3 = vSym.lc3;  lc4 = vSym.lc4;  la4 = vSym.la4;  la2 = vSym.la2;
    K_AFz = vSym.K_AFz;  K_AFx = vSym.K_AFx;  K_AMy = vSym.K_AMy;
    K_BFz = vSym.K_BFz;  K_BFx = vSym.K_BFx;  K_BMy = vSym.K_BMy;
    K_CFz = vSym.K_CFz;  K_CFx = vSym.K_CFx;  K_CMy = vSym.K_CMy;
    K_AFzN3 = vSym.K_AFzN3;  K_AFxN3 = vSym.K_AFxN3;  K_AMyN3 = vSym.K_AMyN3;
    K_BFzN3 = vSym.K_BFzN3;  K_BFxN3 = vSym.K_BFxN3;  K_BMyN3 = vSym.K_BMyN3;
    D_AFz = vSym.D_AFz;  D_AFx = vSym.D_AFx;  D_AMy = vSym.D_AMy;
    D_BFz = vSym.D_BFz;  D_BFx = vSym.D_BFx;  D_BMy = vSym.D_BMy;
    D_CFz = vSym.D_CFz;  D_CFx = vSym.D_CFx;  D_CMy = vSym.D_CMy;
    m1 = vSym.m1;  m2 = vSym.m2;  m3 = vSym.m3;  m4 = vSym.m4;  g = vSym.g;
    I_G1z = vSym.I_G1z;  I_G2z = vSym.I_G2z;  I_G3z = vSym.I_G3z;  I_G4z = vSym.I_G4z;

    % This is for nonlinear interaction spring-damper
    f_intBi = D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4))))) - K_BFxN3*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))))^3;
    f_intBj = - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1))))) - K_BFzN3*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))))^3;
    tau_intBz = K_BMy*(r_q4 + r_q5 - h_q2) + K_BMyN3*(r_q4 + r_q5 - h_q2)^3 + D_BMy*(r_dq4 + r_dq5 - h_dq2);
    f_intAi = K_AFx * r_d2 + K_AFxN3 * r_d2^3 +  D_AFx * r_dd2;
    f_intAj = - K_AFz * r_d3 + K_AFzN3 * r_d3^3 - D_AFz * r_dd3;
    tau_intAz  = K_AMy*r_q4 + K_AMyN3*r_q4^3 + D_AMy*r_dq4;
end

function [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5] = calc_equation_of_motion(vSym, vVar)
    % Assign values to symbols
    l1 = vSym.l1;  l2 = vSym.l2;  l3 = vSym.l3;  l4 = vSym.l4;
    lc1 = vSym.lc1;  la1 = vSym.la1;  lb1 = vSym.lb1;  lc2 = vSym.lc2;
    lc3 = vSym.lc3;  lc4 = vSym.lc4;  la4 = vSym.la4;  la2 = vSym.la2;
    K_AFz = vSym.K_AFz;  K_AFx = vSym.K_AFx;  K_AMy = vSym.K_AMy;
    K_BFz = vSym.K_BFz;  K_BFx = vSym.K_BFx;  K_BMy = vSym.K_BMy;
    K_CFz = vSym.K_CFz;  K_CFx = vSym.K_CFx;  K_CMy = vSym.K_CMy;

    K_AFzN3 = vSym.K_AFzN3;  K_AFxN3 = vSym.K_AFxN3;  K_AMyN3 = vSym.K_AMyN3;
    K_BFzN3 = vSym.K_BFzN3;  K_BFxN3 = vSym.K_BFxN3;  K_BMyN3 = vSym.K_BMyN3;

    D_AFz = vSym.D_AFz;  D_AFx = vSym.D_AFx;  D_AMy = vSym.D_AMy;
    D_BFz = vSym.D_BFz;  D_BFx = vSym.D_BFx;  D_BMy = vSym.D_BMy;
    D_CFz = vSym.D_CFz;  D_CFx = vSym.D_CFx;  D_CMy = vSym.D_CMy;
    m1 = vSym.m1;  m2 = vSym.m2;  m3 = vSym.m3;  m4 = vSym.m4;  g = vSym.g;
    I_G1z = vSym.I_G1z;  I_G2z = vSym.I_G2z;  I_G3z = vSym.I_G3z;  I_G4z = vSym.I_G4z;

    q1 = vVar.q1;  h_q2 = vVar.h_q2;  dq1 = vVar.dq1;  h_dq2 = vVar.h_dq2;
    r_d2 = vVar.r_d2;  r_dd2 = vVar.r_dd2;

    r_d3 = vVar.r_d3;  r_dd3 = vVar.r_dd3;
    r_q4 = vVar.r_q4;  r_dq4 = vVar.r_dq4;
    r_q5 = vVar.r_q5;  r_dq5 = vVar.r_dq5;
    tau1 = vVar.tau1; tau2 = vVar.tau2; tau3 = vVar.tau3; tau4 = vVar.tau4; 

    % nonlienar
    f_intBi = D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4))))) - K_BFxN3*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))))^3;
    f_intBj = - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1))))) - K_BFzN3*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))))^3;

    % linear
    %f_intBi = D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))));
    %f_intBj = - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))));

    f_intCi = D_CFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_CFx*((l2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (l2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + l4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + l4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))));
    f_intCj = - D_CFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_CFz*((l2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (l2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + l4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + l4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))));

    m4gi = m4*g*cos(q1+r_q4+r_q5);
    m4gj = -m4*g*sin(q1+r_q4+r_q5);
    tau_intBz = K_BMy*(r_q4 + r_q5 - h_q2) + K_BMyN3*(r_q4 + r_q5 - h_q2)^3 + D_BMy*(r_dq4 + r_dq5 - h_dq2);
    tau_intCz = D_CMy*(r_dq4 - h_dq2 + r_dq5) + K_CMy*(r_q4 - h_q2 + r_q5);

    f_intAi = K_AFx * r_d2 + K_AFxN3 * r_d2^3 +  D_AFx * r_dd2;
    f_intAj = - K_AFz * r_d3 + K_AFzN3 * r_d3^3 - D_AFz * r_dd3;

    %f_intAi  =K_AFx * r_d2 + D_AFx * r_dd2;
    %f_intAj = - K_AFz * r_d3 - D_AFz * r_dd3;

    m3gi  = m3*g*cos(q1+r_q4);
    m3gj  = -m3*g*sin(q1+r_q4);
    tau_intAz  = K_AMy*r_q4 + K_AMyN3*r_q4^3 + D_AMy*r_dq4;
    %m2gi  = subs(m2*g*cos(q1+h_q2), symbols, [symbol_vals, variable_vals]);
    m2gj = -m2*g*sin(q1+h_q2);
    %m1gi  = subs(m1*g*cos(q1), symbols, [symbol_vals, variable_vals]);
    %m1gj   = subs(-m1*g*sin(q1), symbols, [symbol_vals, variable_vals]);

    a1 = - tau_intBz - tau_intCz - lc2*(f_intBj + f_intCj + m2gj); 
    a2 = (m2*lc2^2 + I_G2z);
    b1 = l3*lc4*m4*sin(r_q5)*r_dq4^2 - l3*la4*m4*sin(r_q5)*r_dq5^2 + tau_intAz + tau_intBz + tau_intCz + f_intAj*lc3 + f_intBj*lc4 + f_intCj*lc4 - lc3*m3gj - lc4*m4gj + f_intBj*l3*cos(r_q5) + f_intCj*l3*cos(r_q5) - l3*m4gj*cos(r_q5) + f_intBi*l3*sin(r_q5) + f_intCi*l3*sin(r_q5) - l3*m4gi*sin(r_q5);
    b2 = - (l3*m4*sin(r_q4) + lc3*m3*sin(r_q4) + lc4*m4*sin(r_q4 + r_q5));
    b3 = - (l3*m4*cos(r_q4) + lc3*m3*cos(r_q4) + lc4*m4*cos(r_q4 + r_q5)); 
    b4 = m4*l3^2 + lc4*m4*cos(r_q5)*l3 + m3*lc3^2 + 2*I_G4z;
    b5 = I_G4z + la4*lc4*m4 + l3*la4*m4*cos(r_q5);
    c1 = l3*lc4*m4*sin(r_q5)*r_dq4^2 + tau_intBz + tau_intCz + f_intBj*lc4 + f_intCj*lc4 - lc4*m4gj;
    c2 = - lc4*m4*sin(r_q4 + r_q5);
    c3 = - lc4*m4*cos(r_q4 + r_q5);
    c4 = I_G4z + l3*lc4*m4*cos(r_q5);
    c5 = (I_G4z + la4*lc4*m4);
    d1 = + f_intAi - m3gi + f_intBi*cos(r_q5) + f_intCi*cos(r_q5) - m4gi*cos(r_q5) - f_intBj*sin(r_q5) - f_intCj*sin(r_q5) + m4gj*sin(r_q5) - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 - la4*m4*r_dq5^2*cos(r_q5);
    d2 = (m3*cos(r_q4) + m4*cos(r_q4));
    d3 = - m3*sin(r_q4) - m4*sin(r_q4);
    d5 = -la4*m4*sin(r_q5);
    e1 = - la4*m4*sin(r_q5)*r_dq5^2 + f_intAj - m3gj + f_intBj*cos(r_q5) + f_intCj*cos(r_q5) - m4gj*cos(r_q5) + f_intBi*sin(r_q5) + f_intCi*sin(r_q5) - m4gi*sin(r_q5);
    e2 = - (m3*sin(r_q4) + m4*sin(r_q4));
    e3 = - (m3*cos(r_q4) + m4*cos(r_q4));
    e4 = l3*m4 + lc3*m3;
    e5 = la4*m4*cos(r_q5);

    h_ddq2 = -(a1 - tau2)/a2;
    r_ddd2 = -(b1*c3*d5*e4 + b1*c4*d3*e5 - b1*c4*d5*e3 - b1*c5*d3*e4 - b3*c1*d5*e4 - b3*c4*d1*e5 + b3*c4*d5*e1 + b3*c5*d1*e4 - b4*c1*d3*e5 + b4*c1*d5*e3 + b4*c3*d1*e5 - b4*c3*d5*e1 - b4*c5*d1*e3 + b4*c5*d3*e1 + b5*c1*d3*e4 - b5*c3*d1*e4 + b5*c4*d1*e3 - b5*c4*d3*e1 + b3*d5*e4*tau4 + b4*d3*e5*tau4 - b4*d5*e3*tau4 - b5*d3*e4*tau4)/(b2*c3*d5*e4 + b2*c4*d3*e5 - b2*c4*d5*e3 - b2*c5*d3*e4 - b3*c2*d5*e4 - b3*c4*d2*e5 + b3*c4*d5*e2 + b3*c5*d2*e4 - b4*c2*d3*e5 + b4*c2*d5*e3 + b4*c3*d2*e5 - b4*c3*d5*e2 - b4*c5*d2*e3 + b4*c5*d3*e2 + b5*c2*d3*e4 - b5*c3*d2*e4 + b5*c4*d2*e3 - b5*c4*d3*e2);
    r_ddd3 = (b1*c2*d5*e4 + b1*c4*d2*e5 - b1*c4*d5*e2 - b1*c5*d2*e4 - b2*c1*d5*e4 - b2*c4*d1*e5 + b2*c4*d5*e1 + b2*c5*d1*e4 - b4*c1*d2*e5 + b4*c1*d5*e2 + b4*c2*d1*e5 - b4*c2*d5*e1 - b4*c5*d1*e2 + b4*c5*d2*e1 + b5*c1*d2*e4 - b5*c2*d1*e4 + b5*c4*d1*e2 - b5*c4*d2*e1 + b2*d5*e4*tau4 + b4*d2*e5*tau4 - b4*d5*e2*tau4 - b5*d2*e4*tau4)/(b2*c3*d5*e4 + b2*c4*d3*e5 - b2*c4*d5*e3 - b2*c5*d3*e4 - b3*c2*d5*e4 - b3*c4*d2*e5 + b3*c4*d5*e2 + b3*c5*d2*e4 - b4*c2*d3*e5 + b4*c2*d5*e3 + b4*c3*d2*e5 - b4*c3*d5*e2 - b4*c5*d2*e3 + b4*c5*d3*e2 + b5*c2*d3*e4 - b5*c3*d2*e4 + b5*c4*d2*e3 - b5*c4*d3*e2);
    r_ddq4 = (b1*c2*d3*e5 - b1*c2*d5*e3 - b1*c3*d2*e5 + b1*c3*d5*e2 + b1*c5*d2*e3 - b1*c5*d3*e2 - b2*c1*d3*e5 + b2*c1*d5*e3 + b2*c3*d1*e5 - b2*c3*d5*e1 - b2*c5*d1*e3 + b2*c5*d3*e1 + b3*c1*d2*e5 - b3*c1*d5*e2 - b3*c2*d1*e5 + b3*c2*d5*e1 + b3*c5*d1*e2 - b3*c5*d2*e1 - b5*c1*d2*e3 + b5*c1*d3*e2 + b5*c2*d1*e3 - b5*c2*d3*e1 - b5*c3*d1*e2 + b5*c3*d2*e1 + b2*d3*e5*tau4 - b2*d5*e3*tau4 - b3*d2*e5*tau4 + b3*d5*e2*tau4 + b5*d2*e3*tau4 - b5*d3*e2*tau4)/(b2*c3*d5*e4 + b2*c4*d3*e5 - b2*c4*d5*e3 - b2*c5*d3*e4 - b3*c2*d5*e4 - b3*c4*d2*e5 + b3*c4*d5*e2 + b3*c5*d2*e4 - b4*c2*d3*e5 + b4*c2*d5*e3 + b4*c3*d2*e5 - b4*c3*d5*e2 - b4*c5*d2*e3 + b4*c5*d3*e2 + b5*c2*d3*e4 - b5*c3*d2*e4 + b5*c4*d2*e3 - b5*c4*d3*e2);
    r_ddq5 = -(b1*c2*d3*e4 - b1*c3*d2*e4 + b1*c4*d2*e3 - b1*c4*d3*e2 - b2*c1*d3*e4 + b2*c3*d1*e4 - b2*c4*d1*e3 + b2*c4*d3*e1 + b3*c1*d2*e4 - b3*c2*d1*e4 + b3*c4*d1*e2 - b3*c4*d2*e1 - b4*c1*d2*e3 + b4*c1*d3*e2 + b4*c2*d1*e3 - b4*c2*d3*e1 - b4*c3*d1*e2 + b4*c3*d2*e1 + b2*d3*e4*tau4 - b3*d2*e4*tau4 + b4*d2*e3*tau4 - b4*d3*e2*tau4)/(b2*c3*d5*e4 + b2*c4*d3*e5 - b2*c4*d5*e3 - b2*c5*d3*e4 - b3*c2*d5*e4 - b3*c4*d2*e5 + b3*c4*d5*e2 + b3*c5*d2*e4 - b4*c2*d3*e5 + b4*c2*d5*e3 + b4*c3*d2*e5 - b4*c3*d5*e2 - b4*c5*d2*e3 + b4*c5*d3*e2 + b5*c2*d3*e4 - b5*c3*d2*e4 + b5*c4*d2*e3 - b5*c4*d3*e2);

end

function T = calc_homotransmtx(a_i_1, alpha_i_1, d_i, theta_i)
    % Calculate the transformation matrix based on DH parameters.
    %
    % Parameters:
    % a_i_1 (float): The parameter a_{i-1}
    % alpha_i_1 (float): The parameter alpha_{i-1} (in radians)
    % d_i (float): The parameter d_i
    % theta_i (float): The parameter theta_i (in radians)
    %
    % Returns:
    % T: The transformation matrix

    % First matrix (D_X(a_{i-1}) * R_X(\alpha_{i-1}))
    cos_alpha = cos(alpha_i_1);
    if ~isa(alpha_i_1, 'sym')
        if abs(cos_alpha) < 1e-10
            cos_alpha = 0;
        end
    end
    matrix_1 = [1, 0, 0, a_i_1;
                0, cos_alpha,      -sin(alpha_i_1), 0;
                0, sin(alpha_i_1), cos_alpha,       0;
                0, 0, 0, 1];
    % Second matrix (D_Z(d_i) * R_Z(theta_i))
    cos_theta = cos(theta_i);
    if ~isa(theta_i, 'sym')
        if abs(cos_theta) < 1e-10
            cos_theta = 0;
        end
    end
    matrix_2 = [cos_theta,    -sin(theta_i), 0, 0;
                sin(theta_i), cos_theta,     0, 0;
                0, 0, 1, d_i;
                0, 0, 0, 1];
    % Matrix multiplication
    T = matrix_1 * matrix_2;
end


function [j1, h_j2, h_j3, r_j2, r_j3, r_j4, r_j5, r_j6, h_cuff_A, h_cuff_B, h_cuff_C, r_cuff_B, T] = calc_joint_position(vSym, q1, h_q2, r_d2, r_d3, r_q4, r_q5)
    H_T01 = calc_homotransmtx(0, 0, 0, -pi/2 + q1);
    H_T12 = calc_homotransmtx(vSym.l1, 0, 0, h_q2);
    H_T23a = calc_homotransmtx(vSym.la2, 0, 0, 0);
    H_T23 = calc_homotransmtx(vSym.l2, 0, 0, 0);

    H_T02 = H_T01 * H_T12;
    H_T03 = H_T02 * H_T23;
    H_T03a = H_T02 * H_T23a;

    R_T01 = calc_homotransmtx(0, 0, 0, q1);
    R_T12 = calc_homotransmtx(0, pi/2, vSym.la1 + r_d2, pi/2);
    R_T23 = calc_homotransmtx(0, -pi/2, r_d3, -pi/2);
    R_T34 = calc_homotransmtx(0, -pi/2, 0, r_q4);
    R_T45 = calc_homotransmtx(vSym.l3, 0, 0, r_q5);
    R_T56a = calc_homotransmtx(vSym.la4, 0, 0, 0);
    R_T56 = calc_homotransmtx(vSym.l4, 0, 0, 0);
    R_T12_zero = calc_homotransmtx(0, pi/2, vSym.la1 + 0, pi/2);

    
    R_T02 = R_T01 * R_T12;
    R_T03 = R_T02 * R_T23;
    R_T04 = R_T03 * R_T34;
    R_T05 = R_T04 * R_T45;
    R_T06a = R_T05 * R_T56a;


    R_T02 = R_T01 * R_T12;
    R_T03 = R_T02 * R_T23;
    R_T04 = R_T03 * R_T34;
    R_T05 = R_T04 * R_T45;
    R_T06 = R_T05 * R_T56;
    R_T06a = R_T05 * R_T56a;
    R_T02_zero = R_T01 * R_T12_zero;

    H_T3a0 = inv_homotransmtx(H_T03a);
    T.H3aR6a_T = H_T3a0*R_T06a;

    j1 = [0; 0; 0];
    h_j2 = H_T02(1:3, 4);
    h_j3 = H_T03(1:3, 4);
    r_j2 = R_T02(1:3, 4);
    r_j3 = R_T03(1:3, 4);
    r_j4 = R_T04(1:3, 4);
    r_j5 = R_T05(1:3, 4);
    r_j6 = R_T06(1:3, 4);
    h_cuff_B = H_T03a(1:3, 4);
    h_cuff_A = R_T02_zero(1:3, 4);
    h_cuff_C = H_T03(1:3, 4);
    r_cuff_B = R_T06a(1:3, 4);
end

function plot_point(points, colors)
    % Plot a list of 3D points.
    %
    % Parameters:
    % points: A cell array where each cell contains a 1x3 array of 3D points.

    numPoints = length(points);

    % Initialize arrays to store the coordinates
    xs = zeros(1, numPoints);
    ys = zeros(1, numPoints);

    % Extract the x, y, and z coordinates
    for i = 1:numPoints
        hold on;
        xs(i) = points{i}(1);
        ys(i) = points{i}(2);
        scatter(xs(i), ys(i), 'LineWidth', 2, 'MarkerEdgeColor', colors(i));

    end

    % Plot each point as a dot
    % scatter3(xs, ys, zs, 'o'); % 'o' is for circle marker

end



function plot_2d_links(points, line_style)
    % Plot a list of 3D points.
    %
    % Parameters:
    % points: A cell array where each cell contains a 1x3 array of 3D points.

    % Number of points
    numPoints = length(points);

    % Initialize arrays to store the coordinates
    xs = zeros(1, numPoints);
    ys = zeros(1, numPoints);

    % Extract the x, y, and z coordinates
    for i = 1:numPoints
        xs(i) = points{i}(1);
        ys(i) = points{i}(2);
    end

    % Plot each point as a dot
    % scatter3(xs, ys, zs, 'o'); % 'o' is for circle marker

    % Hold on to add lines to the same plot
    hold on;

    % Connect the points with lines
    plot(xs, ys, 'LineWidth', 2, 'LineStyle', line_style);

    % Setting labels for axes
    xlabel('X Axis');
    ylabel('Y Axis');

    % Show the plot
    grid on; % Optional: Adds a grid to the plot
end

function plot_3d_points(points)
    % Plot a list of 3D points.
    %
    % Parameters:
    % points: A cell array where each cell contains a 1x3 array of 3D points.

    % Number of points
    numPoints = length(points);

    % Initialize arrays to store the coordinates
    xs = zeros(1, numPoints);
    ys = zeros(1, numPoints);
    zs = zeros(1, numPoints);

    % Extract the x, y, and z coordinates
    for i = 1:numPoints
        xs(i) = points{i}(1);
        ys(i) = points{i}(2);
        zs(i) = points{i}(3);
    end

    % Plot each point as a dot
    % scatter3(xs, ys, zs, 'o'); % 'o' is for circle marker

    % Hold on to add lines to the same plot
    hold on;

    % Connect the points with lines
    plot3(xs, ys, zs, 'LineWidth', 2);

    % Setting labels for axes
    xlabel('X Axis');
    ylabel('Y Axis');
    zlabel('Z Axis');

    % Show the plot
    grid on; % Optional: Adds a grid to the plot
end

function [inv_T] = inv_homotransmtx(T)
    R = T(1:3, 1:3);
    d = T(1:3, 4);
    inv_T = [R.', -R.'*d;
             0,0,0,1];
end
