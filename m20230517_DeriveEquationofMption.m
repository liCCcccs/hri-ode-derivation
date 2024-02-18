
syms l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4
syms q_1 q_2 q_3 q_4
syms d_x d_y

L1 = [l_a1 - l_1; 0; 0; 1];
L2 = [l_a2 - l_2; 0; 0; 1];
L3 = [l_a3 - l_3; 0; 0; 1];
L4 = [l_a4 - l_4; 0; 0; 1];
T10 = [cos(q_1), sin(q_1), 0, cos(q_1)*l_1;
       -sin(q_1), cos(q_1), 0, sin(q_1)*l_1;
       0, 0, 1, 0;
       0, 0, 0, 1];
T03 = [cos(q_3), -sin(q_3), 0, d_x;
       sin(q_3), cos(q_3), 0, d_y;
       0, 0, 1, 0;
       0, 0, 0, 1];
T21 = [cos(q_2), sin(q_2), 0, 0;
       -sin(q_2), cos(q_2), 0, 0;
       0, 0, 1, 0;
       0, 0, 0, 1];
T34 = [cos(q_4), -sin(q_2), 0, 0;
       sin(q_4), cos(q_4), 0, 0;
       0, 0, 1, 0;
       0, 0, 0, 1];

Delta_d_13 = L1 - T10*T03*L3;
simplify(Delta_d_13);

Delta_d_24 = L2 - T21*T10*T03*T34*L4;
simplify(Delta_d_24);

%% State space formulation
vParam.l_a1 = 0.5; vParam.l_1 = 1; vParam.l_c1 = 0.5;
vParam.l_a2 = 0.5; vParam.l_2 = 1; vParam.l_c2 = 0.5;
vParam.l_a3 = 0.5; vParam.l_3 = 1; vParam.l_c3 = 0.5;
vParam.l_a4 = 0.5; vParam.l_4 = 1; vParam.l_c4 = 0.5;
vParam.I_G1z = 1; vParam.I_G2z = 1; vParam.I_G3z = 1; vParam.I_G4z = 1;
vParam.K_t1x = 1; vParam.K_t1y = 1; vParam.K_t2x = 1; vParam.K_t2y = 1;
vParam.K_r1 = 1; vParam.K_r2 = 1;
vParam.m_1 = 1; vParam.m_2 = 1; vParam.m_3 = 1; vParam.m_4 = 1;
vParam.g = 10;
vParam.d_x = 0; vParam.d_y = 0;
vParam.ddD_x = 0; vParam.ddD_y = 0;



syms tau_1 tau_2 tau_3 tau_4
syms F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z
syms ddq_1 ddq_2 ddq_3 ddq_4 dq_1 dq_2 dq_3 dq_4
syms ddD_x ddD_y dD_x dD_y d_x d_y
syms l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4
syms l_c1 l_c2 l_c3 l_c4
syms q_1 q_2 q_3 q_4
syms d_x d_y
syms I_G1z I_G2z I_G3z I_G4z
syms K_t1x K_t1y K_t2x K_t2y K_r1 K_r2
syms m_1 m_2 m_3 m_4
syms g
q_1 = 0; dq_1 = 0; ddq_1 = 0;

eqns = [%tau_1 == tau_2 + m_1*l_c1^2*ddq_1 + m_1*l_c1*g*cos(q_1) + m_2*l_1*g*cos(q_1) - m_2*l_c2*l_1*sin(q_2)*(dq_1+dq_2)^2 + m_2*l_1^2*ddq_1 + m_2*l_c2*l_1*cos(q_2)*(ddq_1+ddq_2 ) - l_1*sin(q_2)*F_int2x - l_1*cos(q_2)*F_int2y + I_G1z*ddq_1;
        tau_2 == I_G2z*(ddq_1 + ddq_2) + m_2*l_c2^2*(ddq_1+ddq_2) + m_2*l_1*l_c2*sin(q_2)*dq_1^2 + m_2*l_1*l_c2*cos(q_2)*ddq_1 + m_2*l_c2*g*cos(q_1 +  q_2) - l_c2*F_int2y - tau_int2z;
        
        tau_3 == tau_4 + m_3*l_c3^2*ddq_3 + m_3*l_c3*ddD_y + m_3*l_c3*g*cos(q_3) + m_4*l_3*g*cos(q_3) - m_4*l_c4*l_3*sin(q_4)*(dq_3 + dq_4)^2 + m_4*l_3^2*(ddq_3 + ddD_y) + l_3*sin(q_4)*F_int2x + l_3*cos(q_4)*F_int2y + I_G3z*ddq_3 + tau_int1z;
        tau_4 == I_G4z*(ddq_3 + ddq_4) + m_4*l_c4^2*(ddq_3 + ddq_4) + m_4*l_3*l_c4*sin(q_4)*(dq_3^2 + ddD_x) + m_4*l_3*l_c4*cos(q_4)*(ddq_3 + ddD_y) + m_4*l_c4*g*cos(q_3 + q_4) + l_c4*F_int2y + tau_int2z;

        0 == F_int1x - sin(q_4)*(F_int2y + l_c4*m_4*(ddq_3 + ddq_4) - g*m_4*cos(q_3 + q_4) + l_3*m_4*cos(q_4)*(ddD_y + ddq_3) + l_3*m_4*sin(q_4)*(dq_3^2 + ddD_x)) + ddD_x*m_3 - cos(q_4)*(l_c4*m_4*(dq_3 + dq_4)^2 - F_int2x + g*m_4*sin(q_3 + q_4) - l_3*m_4*sin(q_4)*(ddD_y + ddq_3) + l_3*m_4*cos(q_4)*(dq_3^2 + ddD_x)) + g*m_3*sin(q_3) - dq_3^2*l_c3*m_3;
        0 ==  F_int1y + ddD_y*m_3 - sin(q_4)*(l_c4*m_4*(dq_3 + dq_4)^2 - F_int2x + g*m_4*sin(q_3 + q_4) - l_3*m_4*sin(q_4)*(ddD_y + ddq_3) + l_3*m_4*cos(q_4)*(dq_3^2 + ddD_x)) + cos(q_4)*(F_int2y + l_c4*m_4*(ddq_3 + ddq_4) - g*m_4*cos(q_3 + q_4) + l_3*m_4*cos(q_4)*(ddD_y + ddq_3) + l_3*m_4*sin(q_4)*(dq_3^2 + ddD_x)) + g*m_3*cos(q_3) + ddq_3*l_c3*m_3;

        F_int1x == K_t1x*(l_a1 - l_1 - d_x*cos(q_1) - l_1*cos(q_1) - d_y*sin(q_1) + cos(q_1 - q_3)*(l_3 - l_a3));
        F_int1y == K_t1y*(d_x*sin(q_1) - d_y*cos(q_1) - sin(q_1 - q_3)*(l_3 - l_a3) - l_1*sin(q_1));
        F_int2x == K_t2x*(l_a2 - l_2 + l_4*cos(q_1 + q_2 - q_3 - q_4) - l_a4*cos(q_1 + q_2 - q_3 - q_4) - d_x*cos(q_1 + q_2) - d_y*sin(q_1 + q_2) - l_1*cos(q_1 - q_2));
        F_int2y == K_t2y*(l_a4*sin(q_1 + q_2 - q_3 - q_4) - l_4*sin(q_1 + q_2 - q_3 - q_4) - d_y*cos(q_1 + q_2) + d_x*sin(q_1 + q_2) - l_1*sin(q_1 - q_2));
        tau_int1z == (q_3 - q_1)*K_r1;
        tau_int2z == (q_3+q_4-q_1-q_2)*K_r2;
        ];

tobeElimitedVars = [F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z];
%tobeElimitedVars = [];
expr = eliminate(eqns, tobeElimitedVars);
ddthetaSolved = solve(expr, [ddq_2, ddq_3, ddq_4, ddD_x, ddD_y]);
%ddthetaSolved = solve(expr, [ddq_2, ddq_3, ddq_4]);


%% solve the dynamic equations
tspan = [0 2];
y0 = [pi/2, 0, pi/2, 0, 0, 0, 0, 0, 0.5, 0.5];  %initial condition
tau = [0, 0, 0, 0];
[t,y] = ode45(@(t,y) runrobot(t,y, ddthetaSolved, tau, vParam), tspan, y0);

%% Animation
l1 = vParam.l_1;
l2 = vParam.l_2;
l3 = vParam.l_3;
l4 = vParam.l_4;
dt = 0.1;
q1 = y(:, 1); q2 = y(:, 2); q3 = y(:, 3); q4 = y(:, 4);
dx = y(:, 9); dy = y(:, 10);

for i = 1:length(y)
    clf
    plot([0, l1*cos(q1(i))],[0, l1*sin(q1(i))],'linewidth',3);  % Human link 1
    hold on
    plot([l1*cos(q1(i)), l1*cos(q1(i))+l2*cos(q1(i)+q2(i))],[l1*sin(q1(i)), l1*sin(q1(i))+l2*sin(q1(i)+q2(i))],'linewidth',1.5)  % Human link 2
 
    plot([dx, dx+l3*cos(q3(i))],[dy, dy+l3*sin(q3(i))],'linewidth',3);  % Robot link 1
    plot([l3*cos(q3(i)), l3*cos(q3(i))+l4*cos(q3(i)+q4(i))],[l3*sin(q3(i)), l3*sin(q3(i))+l4*sin(q3(i)+q4(i))],'linewidth',1.5)  % Robot link 2

    axis([-2,2,-2,2])
    pause(dt)
end


%% Calculate differetial equations
function dydt = runrobot(t, y, ddthetaSolved, tau, vParam)
    % y(1)=q(1), y(2)=q(2), y(3)=q(3), y(4)=q(4), 
    % y(5)=qdot(1), y(6)=qdot(2), y(7)=qdot(3), y(8)=qdot(4)

    % hard limits on joint2 (There is no limits on joints)
%     if y(2) < 0.5
%         y(2) = 0.5;
%         y(4) = 0;
%     elseif y(2) > 1
%         y(2) = 1;
%         y(4) = 0;
%     end
    % x represents the state that contains q and qdot for the robot 
    % ie. for both joints
    syms m_1 m_2 m_3 m_4 
    syms I_G1z I_G2z I_G3z I_G4z 
    syms l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 
    syms l_c1 l_c2 l_c3 l_c4
    syms tau_1  tau_2 tau_3 tau_4
    syms K_t1x K_t1y K_t2x K_t2y K_r1 K_r2
    syms g 
    syms q_1 q_2 q_3 q_4
    syms dq_1 dq_2 dq_3 dq_4
    syms F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z
    syms d_x d_y ddD_x ddD_y

    v_theta1 = y(1); v_theta2 = y(2); v_theta3 = y(3); v_theta4 = y(4);
    v_dtheta1 = y(5); v_dtheta2 = y(6); v_dtheta3 = y(7); v_dtheta4 = y(8);
    v_d_x = y(9); v_d_y = y(10); v_dD_x = y(9); 

    v_tau_1 = tau(1); v_tau_2 = tau(2); v_tau_3 = tau(3); v_tau_4 = tau(4);

        v_F_int1x = subs(K_t1x*(l_a1 - l_1 - d_x*cos(q_1) - l_1*cos(q_1) - d_y*sin(q_1) + cos(q_1 - q_3)*(l_3 - l_a3)), {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
        v_F_int1y = subs(K_t1y*(d_x*sin(q_1) - d_y*cos(q_1) - sin(q_1 - q_3)*(l_3 - l_a3) - l_1*sin(q_1)), {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
        v_F_int2x = subs(K_t2x*(l_a2 - l_2 + l_4*cos(q_1 + q_2 - q_3 - q_4) - l_a4*cos(q_1 + q_2 - q_3 - q_4) - d_x*cos(q_1 + q_2) - d_y*sin(q_1 + q_2) - l_1*cos(q_1 - q_2)), {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
        v_F_int2y = subs(K_t2y*(l_a4*sin(q_1 + q_2 - q_3 - q_4) - l_4*sin(q_1 + q_2 - q_3 - q_4) - d_y*cos(q_1 + q_2) + d_x*sin(q_1 + q_2) - l_1*sin(q_1 - q_2)), {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
        v_tau_int1z = subs((q_3 - q_1)*K_r1, {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
        v_tau_int2z = subs((q_3+q_4-q_1-q_2)*K_r2, {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
  
    ddtheta1 = subs(ddthetaSolved.ddq_1, {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
    ddtheta2 = subs(ddthetaSolved.ddq_2, {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
    ddtheta3 = subs(ddthetaSolved.ddq_3, {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
    ddtheta4 = subs(ddthetaSolved.ddq_4, {m_1 m_2 m_3 m_4 I_G1z I_G2z I_G3z I_G4z l_a1 l_1 l_a2 l_2 l_a3 l_3 l_a4 l_4 l_c1 l_c2 l_c3 l_c4 K_t1x K_t1y K_t2x K_t2y K_r1 K_r2 g q_1 q_2 q_3 q_4 dq_1 dq_2 dq_3 dq_4 tau_1  tau_2 tau_3 tau_4 d_x d_y ddD_x ddD_y}, ...
        {vParam.m_1 vParam.m_2 vParam.m_3 vParam.m_4 vParam.I_G1z vParam.I_G2z vParam.I_G3z vParam.I_G4z vParam.l_a1 vParam.l_1 vParam.l_a2 vParam.l_2 vParam.l_a3 vParam.l_3 vParam.l_a4 vParam.l_4 vParam.l_c1 vParam.l_c2 vParam.l_c3 vParam.l_c4 vParam.K_t1x vParam.K_t1y vParam.K_t2x vParam.K_t2y vParam.K_r1 vParam.K_r2 vParam.g, v_theta1 v_theta2 v_theta3 v_theta4 v_dtheta1 v_dtheta2 v_dtheta3 v_dtheta4 v_tau_1  v_tau_2 v_tau_3 v_tau_4 vParam.d_x vParam.d_y vParam.ddD_x  vParam.ddD_y});
    
    ddtheta1 = subs(ddtheta1, {F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z}, {v_F_int1x v_F_int1y v_F_int2x v_F_int2y v_tau_int1z v_tau_int2z});
    ddtheta2 = subs(ddtheta2, {F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z}, {v_F_int1x v_F_int1y v_F_int2x v_F_int2y v_tau_int1z v_tau_int2z});
    ddtheta3 = subs(ddtheta3, {F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z}, {v_F_int1x v_F_int1y v_F_int2x v_F_int2y v_tau_int1z v_tau_int2z});
    ddtheta4 = subs(ddtheta4, {F_int1x F_int1y F_int2x F_int2y tau_int1z tau_int2z}, {v_F_int1x v_F_int1y v_F_int2x v_F_int2y v_tau_int1z v_tau_int2z});

    

    v_ddtheta1 = double(ddtheta1);
    v_ddtheta2 = double(ddtheta2);
    v_ddtheta3 = double(ddtheta3);
    v_ddtheta4 = double(ddtheta4);

    dydt = zeros(8,1);
    dydt(1) = v_dtheta1; 
    dydt(2) = v_dtheta2;
    dydt(3) = v_dtheta3; 
    dydt(4) = v_dtheta4;
    dydt(5) = v_ddtheta1; 
    dydt(6) = v_ddtheta2;
    dydt(7) = v_ddtheta3;
    dydt(8) = v_ddtheta4;
end





