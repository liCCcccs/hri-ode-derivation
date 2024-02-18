vParam.m1 = 2; %(kg)
vParam.m2 = 1; %(kg)
vParam.IG2y = 0.5; 
vParam.IG1y = 0.3;
vParam.l1 = 1;
vParam.l2 = 0.6;
vParam.g = 9.8;

syms m1 l1 m2 l2 g
syms IG1x IG1y IG1z IG2x IG2y IG2z
syms theta1 dtheta1 ddtheta1 theta2 dtheta2 ddtheta2
syms RA1x RA1y RA1z RB1x RB1y RB1z
syms RA2x RA2y RA2z
syms MA1x MA1y MA1z MB1x MB1y MB1z
syms MA2x MA2y MA2z

eqns = [-m1*l1/2*ddtheta1 == RA1x + RB1x + m1*g*sin(theta1); 
        0 == RA1y + RB1y;
        m1*l1/2*dtheta1^2 == RA1z + RB1z - m1*g*cos(theta1);
        0 == MA1x + MB1x - l1/2*RA1y + l1/2*RB1y;
        IG1y*ddtheta1 == MA1y + MB1y + l1/2*RA1x - l1/2*RB1x;
        0 == MA1z + MB1z;
        -m2*l2/2*(ddtheta1+ddtheta2) == RA2x + m2*g*sin(theta1+theta2);
        0 == RA2y;
        m2*l2/2*(dtheta1+dtheta2)^2 == RA2z - m2*g*cos(theta1+theta2);
        0 == MA2x - l2/2*RA2y;
        IG2y*(ddtheta1+ddtheta2) == MA2y + l2/2*RA2x;
        0 == MA2z;
        RB1x == -cos(theta2)*RA2x - sin(theta2)*RA2z;
        RB1y == -RA2y;
        RB1z == sin(theta2)*RA2x - cos(theta2)*RA2z;
        MB1x == -cos(theta2)*MA2x - sin(theta2)*MA2z;
        MB1y == -MA2y;
        MB1z == sin(theta2)*MA2x - cos(theta2)*MA2z
        ];

tobeElimitedVars = [RA1x RA1y RA1z RB1x RB1y RB1z RA2x RA2y RA2z MA1x MA1z MB1x MB1y MB1z MA2x MA2z];
expr = eliminate(eqns, tobeElimitedVars);

ddthetaSolved = solve(expr, [ddtheta1, ddtheta2]);
tauSolved = solve(expr, [MA1y, MA2y]);

%% solve the dynamic equations
tspan = [0 5];
y0 = [pi/2, 0, 0, 0];  %initial condition
tau = [0, 0];
[t,y] = ode45(@(t,y) runrobot(t,y, ddthetaSolved, tau, vParam), tspan, y0);

%% Animation
l1 = 1;
l2 = 0.6;
dt = 0.1;
y1 = y(:, 1) - pi/2; y2 = y(:, 2);
for i = 1:length(y)
    clf
    plot([0, l1*cos(y1(i))],[0, l1*sin(y1(i))],'linewidth',3);
    hold on
    plot([l1*cos(y1(i)), l1*cos(y1(i))+l2*cos(y1(i)+y2(i))],[l1*sin(y1(i)), l1*sin(y1(i))+l2*sin(y1(i)+y2(i))],'linewidth',1.5)
    axis([-2,2,-2,2])
    pause(dt)
end


%% Calculate differetial equations
function dydt = runrobot(t, y, ddthetaSolved, tau, vParam)
    %y(1)=q(1), y(2)=q(2), y(3)=qdot(1), y(4)=qdot(2), 
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
    syms m1 m2 IG1y IG2y l1 l2 g MA1y MA2y theta1 theta2 dtheta1 dtheta2

    v_theta1 = y(1);
    v_theta2 = y(2);
    v_dtheta1 = y(3);
    v_dtheta2 = y(4);
    v_MA1y = tau(1);
    v_MA2y = tau(2);
  
    %ddtheta1 = (2*(8*IG2y*MA1y - 8*IG2y*MA2y + 2*MA1y*l2^2*m2 - 2*MA2y*l2^2*m2 + dtheta1^2*l1*l2^3*m2^2*sin(theta2) + dtheta2^2*l1*l2^3*m2^2*sin(theta2) - 4*MA2y*l1*l2*m2*cos(theta2) - 4*IG2y*g*l1*m1*sin(theta1) + 2*dtheta1*dtheta2*l1*l2^3*m2^2*sin(theta2) + 2*g*l1*l2^2*m2^2*cos(theta1 + theta2)*sin(theta2) + 8*IG2y*g*l1*m2*cos(theta1 + theta2)*sin(theta2) - 8*IG2y*g*l1*m2*sin(theta1 + theta2)*cos(theta2) + 4*IG2y*dtheta1^2*l1*l2*m2*sin(theta2) + 4*IG2y*dtheta2^2*l1*l2*m2*sin(theta2) - g*l1*l2^2*m1*m2*sin(theta1) + 8*IG2y*dtheta1*dtheta2*l1*l2*m2*sin(theta2)))/(m1*m2*l1^2*l2^2 + 4*IG2y*m1*l1^2 + 4*IG1y*m2*l2^2 + 16*IG1y*IG2y);
    %ddtheta2 = -(2*(8*IG2y*MA1y - 8*IG1y*MA2y - 8*IG2y*MA2y - 2*MA2y*l1^2*m1 + 2*MA1y*l2^2*m2 - 2*MA2y*l2^2*m2 + dtheta1^2*l1*l2^3*m2^2*sin(theta2) + dtheta2^2*l1*l2^3*m2^2*sin(theta2) + 4*IG1y*g*l2*m2*sin(theta1 + theta2) - 4*MA2y*l1*l2*m2*cos(theta2) - 4*IG2y*g*l1*m1*sin(theta1) + 2*dtheta1*dtheta2*l1*l2^3*m2^2*sin(theta2) + 2*g*l1*l2^2*m2^2*cos(theta1 + theta2)*sin(theta2) + 8*IG2y*g*l1*m2*cos(theta1 + theta2)*sin(theta2) - 8*IG2y*g*l1*m2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*l2*m1*m2*sin(theta1 + theta2) + 4*IG2y*dtheta1^2*l1*l2*m2*sin(theta2) + 4*IG2y*dtheta2^2*l1*l2*m2*sin(theta2) - g*l1*l2^2*m1*m2*sin(theta1) + 8*IG2y*dtheta1*dtheta2*l1*l2*m2*sin(theta2)))/(m1*m2*l1^2*l2^2 + 4*IG2y*m1*l1^2 + 4*IG1y*m2*l2^2 + 16*IG1y*IG2y);
    ddtheta1 = subs(ddthetaSolved.ddtheta1, {m1, m2, IG1y, IG2y, l1, l2, g, MA1y, MA2y, theta1, theta2, dtheta1, dtheta2}, {vParam.m1, vParam.m2, vParam.IG1y, vParam.IG2y, vParam.l1, vParam.l2, vParam.g, v_MA1y, v_MA2y, v_theta1, v_theta2, v_dtheta1, v_dtheta2});
    ddtheta2 = subs(ddthetaSolved.ddtheta2, {m1, m2, IG1y, IG2y, l1, l2, g, MA1y, MA2y, theta1, theta2, dtheta1, dtheta2}, {vParam.m1, vParam.m2, vParam.IG1y, vParam.IG2y, vParam.l1, vParam.l2, vParam.g, v_MA1y, v_MA2y, v_theta1, v_theta2, v_dtheta1, v_dtheta2});
    v_ddtheta1 = double(ddtheta1);
    v_ddtheta2 = double(ddtheta2);

    dydt = zeros(4,1);
    dydt(1) = v_dtheta1; 
    dydt(2) = v_dtheta2;
    dydt(3) = v_ddtheta1;
    dydt(4) = v_ddtheta2;
end



