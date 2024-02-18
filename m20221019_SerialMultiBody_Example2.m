vParam.mR1 = 2; %(kg)
vParam.mH1 = 1; %(kg)
vParam.IRG1y = 0.5; 
vParam.IHG1y = 0.3;
vParam.lR1 = 1;
vParam.lH1 = 1;
vParam.lRc1 = 0.8;
vParam.lHc1 = 0.8;
vParam.g = 9.8;
vParam.KFx = 1;
vParam.KFy = 1;
vParam.KFz = 1;
vParam.KMx = 1;
vParam.KMy = 1;
vParam.KMz = 1;

syms mR1 lR1 mH1 lH1 g lRc1 lHc1
syms IRG1x IRG1y IRG1z IHG1x IHG1y IHG1z
syms thetaR1 dthetaR1 ddthetaR1 thetaH1 dthetaH1 ddthetaH1
syms RRA1x RRA1y RRA1z FRc1x FRc1y FRc1z
syms RHA1x RHA1y RHA1z FHc1x FHc1y FHc1z
syms MRA1x MRA1y MRA1z MRc1x MRc1y MRc1z
syms MHA1x MHA1y MHA1z MHc1x MHc1y MHc1z
syms KFx KFy KFz KMx KMy KMz

eqns = [-mR1*lR1/2*ddthetaR1 == RRA1x + FRc1x + mR1*g*sin(thetaR1); 
        0 == RRA1y + FRc1y;
        mR1*lR1/2*dthetaR1^2 == RRA1z + FRc1z - mR1*g*cos(thetaR1);
        0 == MRA1x + MRc1x - lR1/2*RRA1y - (lR1/2-lRc1)*FRc1y;
        IRG1y*ddthetaR1 == MRA1y + MRc1y + lR1/2*RRA1x + (lR1/2-lRc1)*FRc1x;
        0 == MRA1z + MRc1z;
        -mH1*lH1/2*(ddthetaH1) == RHA1x + mH1*g*sin(thetaH1) + FHc1x;
        0 == RHA1y + FHc1y;
        mH1*lH1/2*(dthetaH1)^2 == RHA1z - mH1*g*cos(thetaH1) + FHc1z;
        0 == MHA1x + MHc1x - lH1/2*RHA1y - (lH1/2-lHc1)*FHc1y;
        IHG1y*(ddthetaH1) == MHA1y + MHc1y + lH1/2*RHA1x + (lH1/2-lHc1)*FHc1x;
        0 == MHA1z + MHc1z;
        FRc1x == -cos(thetaH1-thetaR1)*FHc1x - sin(thetaH1-thetaR1)*FHc1z;
        FRc1y == -FHc1y;
        FRc1z == sin(thetaH1-thetaR1)*FHc1x - cos(thetaH1-thetaR1)*FHc1z;
        MRc1x == -cos(thetaH1-thetaR1)*MHc1x - sin(thetaH1-thetaR1)*MHc1z;
        MRc1y == -MHc1y;
        MRc1z == sin(thetaH1-thetaR1)*MHc1x - cos(thetaH1-thetaR1)*MHc1z;
        FRc1x == -KFx*sin(thetaH1-thetaR1)*lRc1;
        FRc1y == 0;
        FRc1z == KFz*(lRc1+cos(thetaH1-thetaR1)*lRc1);
        MRc1x == KMx*atan2(0, cos(thetaH1-thetaR1));
        MRc1y == KMy*atan2(sin(thetaH1-thetaR1), cos(thetaH1-thetaR1));
        MRc1z == KMz*atan2(0, cos(thetaH1-thetaR1))
        ];

tobeElimitedVars = [RRA1x RRA1y RRA1z FRc1x FRc1y FRc1z RHA1x RHA1y RHA1z FHc1x FHc1y FHc1z MRA1x MRA1z MRc1x MRc1y MRc1z MHA1x MHA1z MHc1x MHc1y MHc1z];
expr = eliminate(eqns, tobeElimitedVars);

ddthetaSolved = solve(expr, [ddthetaR1, ddthetaH1]);
tauSolved = solve(expr, [MRA1y, MHA1y]);

%% solve the dynamic equations
tspan = [0 5];
y0 = [pi/2, pi/2, 0, 0];  %initial condition
tau = [0, 0];
[t,y] = ode45(@(t,y) runrobot(t,y, ddthetaSolved, tau, vParam), tspan, y0);

%% Animation
l1 = 1;
l2 = 0.6;
dt = 0.1;
y1 = y(:, 1) - pi/2; y2 = y(:, 2) - pi/2;
for i = 1:length(y)
    clf
    plot([0, l1*cos(y1(i))],[0, l1*sin(y1(i))],'linewidth',3);
    hold on
    plot([0, l1*cos(y2(i))],[0, l1*sin(y2(i))],'linewidth',3);
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
    syms mR1 mH1 IRG1y IHG1y lR1 lH1 lRc1 lHc1 g KFx KFy KFz KMx KMy KMz MRA1y MHA1y thetaR1 thetaH1 dthetaR1 dthetaH1
    
    v_thetaR1 = y(1);
    v_thetaH1 = y(2);
    v_dthetaR1 = y(3);
    v_dthetaH1 = y(4);
    v_MRA1y = tau(1);
    v_MHA1y = tau(2);
  
    ddthetaR1 = subs(ddthetaSolved.ddthetaR1, ...
        {mR1, mH1, IRG1y, IHG1y, lR1, lH1, lRc1, lHc1, g, KFx, KFy, KFz, KMx, KMy, KMz, MRA1y, MHA1y, thetaR1, thetaH1, dthetaR1, dthetaH1}, ...
        {vParam.mR1, vParam.mH1, vParam.IRG1y, vParam.IHG1y, vParam.lR1, vParam.lH1, vParam.lRc1, vParam.lHc1, vParam.g, vParam.KFx, vParam.KFy, vParam.KFz, vParam.KMx, vParam.KMy, vParam.KMz, v_MRA1y, v_MHA1y, v_thetaR1, v_thetaH1, v_dthetaR1, v_dthetaH1});
    ddthetaH1 = subs(ddthetaSolved.ddthetaH1, ...
        {mR1, mH1, IRG1y, IHG1y, lR1, lH1, lRc1, lHc1, g, KFx, KFy, KFz, KMx, KMy, KMz, MRA1y, MHA1y, thetaR1, thetaH1, dthetaR1, dthetaH1}, ...
        {vParam.mR1, vParam.mH1, vParam.IRG1y, vParam.IHG1y, vParam.lR1, vParam.lH1, vParam.lRc1, vParam.lHc1, vParam.g, vParam.KFx, vParam.KFy, vParam.KFz, vParam.KMx, vParam.KMy, vParam.KMz, v_MRA1y, v_MHA1y, v_thetaR1, v_thetaH1, v_dthetaR1, v_dthetaH1});
    v_ddthetaR1 = double(ddthetaR1);
    v_ddthetaH1 = double(ddthetaH1);

    dydt = zeros(4,1);
    dydt(1) = v_dthetaR1; 
    dydt(2) = v_dthetaH1;
    dydt(3) = v_ddthetaR1;
    dydt(4) = v_ddthetaH1;
end

















