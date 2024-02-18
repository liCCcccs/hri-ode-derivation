syms l1 l2 l3 l4
syms lc1 la1 lb1 lc2 lc3 lc4 la4 la2
syms q1 h_q2 dq1 h_dq2 ddq1 h_ddq2 
syms r_d2 r_dd2 r_ddd2
syms r_d3 r_dd3 r_ddd3
syms r_q4 r_dq4 r_ddq4
syms r_q5 r_dq5 r_ddq5
syms K_AFz K_AFx K_AMy
syms K_BFz K_BFx K_BMy
syms D_AFz D_AFx D_AMy
syms D_BFz D_BFx D_BMy
syms m1 m2 m3 m4 g
syms I_G1z I_G2z I_G3z I_G4z
syms tau1 tau2 tau3 tau4

syms shq2 srd2 srd3 srq4 srq5
syms chq2 crd2 crd3 crq4 crq5

q1 = pi/2; dq1 = 0; ddq1 = 0;

eqns = [
        tau2 == D_BMy*h_dq2 + I_G2z*ddq1 + I_G2z*h_ddq2 + K_BMy*h_q2 - D_BMy*r_dq4 - D_BMy*r_dq5 - K_BMy*r_q4 - K_BMy*r_q5 + h_ddq2*lc2^2*m2 + ddq1*lc2^2*m2 - K_BFz*la4*lc2*sin(r_q4 - h_q2 + r_q5) + g*lc2*m2*sin(h_q2 + q1) + D_BFz*lc2*r_dd2*chq2 + D_BFz*lc2*r_dd3*chq2 + K_BFz*lc2*r_d2*chq2 - K_BFz*l1*lc2*shq2 + K_BFz*la1*lc2*shq2 + K_BFz*lc2*r_d3*shq2 + K_BFz*l3*lc2*sin(h_q2 - r_q4) + dq1^2*l1*lc2*m2*shq2 + D_BFz*h_dq2*l1*lc2*chq2 - D_BFz*la1*lc2*r_dq4*chq2 - D_BFz*la1*lc2*r_dq5*chq2 - D_BFz*lc2*r_d3*r_dq4*chq2 - D_BFz*lc2*r_d3*r_dq5*chq2 + D_BFz*lc2*r_d2*r_dq4*shq2 + D_BFz*lc2*r_d2*r_dq5*shq2 + ddq1*l1*lc2*m2*chq2 - D_BFz*l3*lc2*r_dq5*cos(h_q2 - r_q4)

        tau3 == 2*I_G4z*ddq1 - D_BMy*h_dq2 - K_BMy*h_q2 + D_AMy*r_dq4 + D_BMy*r_dq4 + D_BMy*r_dq5 + 2*I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_AMy*r_q4 + K_BMy*r_q4 + K_BMy*r_q5 + l3^2*m4*r_ddq4 + lc3^2*m3*r_ddq4 - D_AFz*lc3*r_dd2 - K_AFz*lc3*r_d2 + la4*lc4*m4*r_ddq5 - D_BFz*lc4*r_dd2*chq2 - D_BFz*lc4*r_dd3*chq2 - K_BFz*lc4*r_d2*chq2 + K_BFz*l1*lc4*shq2 - K_BFz*la1*lc4*shq2 - K_BFz*lc4*r_d3*shq2 - K_BFx*l3*la2*srq5 - l3*m4*r_ddd2*crq4 - lc3*m3*r_ddd2*crq4 - l3*m4*r_ddd3*srq4 - lc3*m3*r_ddd3*srq4 + K_BFx*l3^2*shq2*srq4*srq5 - D_BFz*l3*r_dd2*chq2*crq5 - D_BFz*l3*r_dd3*chq2*crq5 - K_BFz*l3*r_d2*chq2*crq5 - K_BFx*l1*l3*chq2*srq5 + K_BFz*l1*l3*crq5*shq2 + K_BFx*l3*la1*chq2*srq5 - K_BFx*l3*la4*chq2*srq4 + K_BFx*l3*la4*crq4*shq2 - K_BFz*l3*la1*crq5*shq2 + K_BFz*l3*lc4*chq2*srq4 - K_BFz*l3*lc4*crq4*shq2 + K_BFx*l3*r_d3*chq2*srq5 - K_BFz*l3*r_d3*crq5*shq2 + dq1^2*l3*la1*m4*crq4 + dq1^2*la1*lc3*m3*crq4 - D_BFx*l3*r_dd2*shq2*srq5 - D_BFx*l3*r_dd3*shq2*srq5 - K_BFx*l3*r_d2*shq2*srq5 - lc4*m4*r_ddd2*crq4*crq5 + g*l3*m4*(0)*srq4 + g*l3*m4*crq4*(-1) + g*lc3*m3*(0)*srq4 + g*lc3*m3*crq4*(-1) - l3*la4*m4*r_dq5^2*srq5 + l3*lc4*m4*r_dq4^2*srq5 - lc4*m4*r_ddd3*crq4*srq5 - lc4*m4*r_ddd3*crq5*srq4 + lc4*m4*r_ddd2*srq4*srq5 - D_BFz*h_dq2*l1*lc4*chq2 + D_BFz*la1*lc4*r_dq4*chq2 + D_BFz*la1*lc4*r_dq5*chq2 + D_BFz*lc4*r_d3*r_dq4*chq2 + D_BFz*lc4*r_d3*r_dq5*chq2 - D_BFz*lc4*r_d2*r_dq4*shq2 - D_BFz*lc4*r_d2*r_dq5*shq2 + l3*la4*m4*r_ddq5*crq5 + l3*lc4*m4*r_ddq4*crq5 + ddq1*l3*la1*m4*srq4 + ddq1*la1*lc3*m3*srq4 + K_BFx*l3^2*chq2*crq4*srq5 + K_BFz*l3^2*chq2*crq5*srq4 - K_BFz*l3^2*crq4*crq5*shq2 + K_BFz*la4*lc4*chq2*crq4*srq5 + K_BFz*la4*lc4*chq2*crq5*srq4 - K_BFz*la4*lc4*crq4*crq5*shq2 + dq1^2*la1*lc4*m4*crq4*crq5 + K_BFz*la4*lc4*shq2*srq4*srq5 + g*lc4*m4*(0)*crq4*srq5 + g*lc4*m4*(0)*crq5*srq4 + g*lc4*m4*crq4*crq5*(-1) - dq1^2*la1*lc4*m4*srq4*srq5 - g*lc4*m4*(-1)*srq4*srq5 + D_BFz*l3^2*r_dq5*chq2*crq4*crq5 + K_BFx*l3*la4*chq2*crq5^2*srq4 - K_BFx*l3*la4*crq4*crq5^2*shq2 + K_BFz*l3*la4*chq2*crq5^2*srq4 - K_BFz*l3*la4*crq4*crq5^2*shq2 - D_BFx*l3^2*r_dq5*chq2*srq4*srq5 + D_BFx*l3^2*r_dq5*crq4*shq2*srq5 + D_BFz*l3^2*r_dq5*crq5*shq2*srq4 - D_BFz*h_dq2*l1*l3*chq2*crq5 + D_BFz*l3*la1*r_dq4*chq2*crq5 + D_BFz*l3*la1*r_dq5*chq2*crq5 + D_BFz*l3*lc4*r_dq5*chq2*crq4 + D_BFz*l3*r_d3*r_dq4*chq2*crq5 + D_BFz*l3*r_d3*r_dq5*chq2*crq5 - D_BFx*h_dq2*l1*l3*shq2*srq5 + D_BFx*l3*r_d2*r_dq4*chq2*srq5 + D_BFx*l3*r_d2*r_dq5*chq2*srq5 - D_BFz*l3*r_d2*r_dq4*crq5*shq2 - D_BFz*l3*r_d2*r_dq5*crq5*shq2 + D_BFx*l3*la1*r_dq4*shq2*srq5 + D_BFx*l3*la1*r_dq5*shq2*srq5 + D_BFz*l3*lc4*r_dq5*shq2*srq4 + D_BFx*l3*r_d3*r_dq4*shq2*srq5 + D_BFx*l3*r_d3*r_dq5*shq2*srq5 + ddq1*la1*lc4*m4*crq4*srq5 + ddq1*la1*lc4*m4*crq5*srq4 + K_BFx*l3*la4*crq5*shq2*srq4*srq5 + K_BFz*l3*la4*crq5*shq2*srq4*srq5 + K_BFx*l3*la4*chq2*crq4*crq5*srq5 + K_BFz*l3*la4*chq2*crq4*crq5*srq5

        tau4 == I_G4z*ddq1 - D_BMy*h_dq2 - K_BMy*h_q2 + D_BMy*r_dq4 + D_BMy*r_dq5 + I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_BMy*r_q4 + K_BMy*r_q5 + K_BFz*la4*lc4*sin(r_q4 - h_q2 + r_q5) + la4*lc4*m4*r_ddq5 - D_BFz*lc4*r_dd2*chq2 - D_BFz*lc4*r_dd3*chq2 - lc4*m4*r_ddd2*cos(r_q4 + r_q5) - K_BFz*lc4*r_d2*chq2 + K_BFz*l1*lc4*shq2 - K_BFz*la1*lc4*shq2 - lc4*m4*r_ddd3*sin(r_q4 + r_q5) - K_BFz*lc4*r_d3*shq2 - K_BFz*l3*lc4*sin(h_q2 - r_q4) + g*lc4*m4*sin(q1 + r_q4 + r_q5) + dq1^2*la1*lc4*m4*cos(r_q4 + r_q5) + l3*lc4*m4*r_dq4^2*srq5 - D_BFz*h_dq2*l1*lc4*chq2 + D_BFz*la1*lc4*r_dq4*chq2 + D_BFz*la1*lc4*r_dq5*chq2 + D_BFz*lc4*r_d3*r_dq4*chq2 + D_BFz*lc4*r_d3*r_dq5*chq2 + ddq1*la1*lc4*m4*sin(r_q4 + r_q5) - D_BFz*lc4*r_d2*r_dq4*shq2 - D_BFz*lc4*r_d2*r_dq5*shq2 + D_BFz*l3*lc4*r_dq5*cos(h_q2 - r_q4) + l3*lc4*m4*r_ddq4*crq5

        0 == D_AFx*r_dd3 + K_AFx*r_d3 - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 + m3*r_ddd3*crq4 + m4*r_ddd3*crq4 - m3*r_ddd2*srq4 - m4*r_ddd2*srq4 - K_BFx*la2*crq5 - ddq1*la1*m3*crq4 - ddq1*la1*m4*crq4 - la4*m4*r_ddq5*srq5 - K_BFx*l1*chq2*crq5 + K_BFx*la1*chq2*crq5 - K_BFz*la4*chq2*crq4 + K_BFx*r_d3*chq2*crq5 - D_BFx*r_dd2*crq5*shq2 - D_BFx*r_dd3*crq5*shq2 + D_BFz*r_dd2*chq2*srq5 + D_BFz*r_dd3*chq2*srq5 - K_BFx*r_d2*crq5*shq2 + K_BFz*r_d2*chq2*srq5 - K_BFz*l1*shq2*srq5 + K_BFz*la1*shq2*srq5 - K_BFz*la4*shq2*srq4 - g*m3*(0)*crq4 - g*m4*(0)*crq4 + K_BFz*r_d3*shq2*srq5 - la4*m4*r_dq5^2*crq5 + dq1^2*la1*m3*srq4 + dq1^2*la1*m4*srq4 + g*m3*(-1)*srq4 + g*m4*(-1)*srq4 + K_BFx*la4*crq5^2*shq2*srq4 + K_BFz*la4*crq5^2*shq2*srq4 - D_BFx*h_dq2*l1*crq5*shq2 + D_BFz*h_dq2*l1*chq2*srq5 + D_BFx*r_d2*r_dq4*chq2*crq5 + D_BFx*r_d2*r_dq5*chq2*crq5 + D_BFx*la1*r_dq4*crq5*shq2 + D_BFx*la1*r_dq5*crq5*shq2 - D_BFz*la1*r_dq4*chq2*srq5 - D_BFz*la1*r_dq5*chq2*srq5 + D_BFx*r_d3*r_dq4*crq5*shq2 + D_BFx*r_d3*r_dq5*crq5*shq2 - D_BFz*r_d3*r_dq4*chq2*srq5 - D_BFz*r_d3*r_dq5*chq2*srq5 + D_BFz*r_d2*r_dq4*shq2*srq5 + D_BFz*r_d2*r_dq5*shq2*srq5 + K_BFx*l3*chq2*crq4*crq5 + K_BFx*l3*crq5*shq2*srq4 - K_BFz*l3*chq2*srq4*srq5 + K_BFz*l3*crq4*shq2*srq5 + K_BFx*la4*chq2*crq4*crq5^2 + K_BFz*la4*chq2*crq4*crq5^2 - D_BFx*l3*r_dq5*chq2*crq5*srq4 + D_BFx*l3*r_dq5*crq4*crq5*shq2 - D_BFz*l3*r_dq5*chq2*crq4*srq5 - D_BFz*l3*r_dq5*shq2*srq4*srq5 - K_BFx*la4*chq2*crq5*srq4*srq5 + K_BFx*la4*crq4*crq5*shq2*srq5 - K_BFz*la4*chq2*crq5*srq4*srq5 + K_BFz*la4*crq4*crq5*shq2*srq5
        0 == l3*m4*r_ddq4 - K_AFz*r_d2 - m3*r_ddd2*crq4 - m4*r_ddd2*crq4 - m3*r_ddd3*srq4 - m4*r_ddd3*srq4 - D_AFz*r_dd2 + lc3*m3*r_ddq4 - K_BFx*la2*srq5 + la4*m4*r_ddq5*crq5 + ddq1*la1*m3*srq4 + ddq1*la1*m4*srq4 - D_BFz*r_dd2*chq2*crq5 - D_BFz*r_dd3*chq2*crq5 - K_BFz*r_d2*chq2*crq5 - K_BFx*l1*chq2*srq5 + K_BFz*l1*crq5*shq2 + K_BFx*la1*chq2*srq5 - K_BFx*la4*chq2*srq4 + K_BFx*la4*crq4*shq2 - K_BFz*la1*crq5*shq2 + K_BFx*r_d3*chq2*srq5 - K_BFz*r_d3*crq5*shq2 + dq1^2*la1*m3*crq4 + dq1^2*la1*m4*crq4 - D_BFx*r_dd2*shq2*srq5 - D_BFx*r_dd3*shq2*srq5 - K_BFx*r_d2*shq2*srq5 + g*m3*(0)*srq4 + g*m3*crq4*(-1) + g*m4*(0)*srq4 + g*m4*crq4*(-1) - la4*m4*r_dq5^2*srq5 - D_BFz*h_dq2*l1*chq2*crq5 + D_BFz*la1*r_dq4*chq2*crq5 + D_BFz*la1*r_dq5*chq2*crq5 + D_BFz*r_d3*r_dq4*chq2*crq5 + D_BFz*r_d3*r_dq5*chq2*crq5 - D_BFx*h_dq2*l1*shq2*srq5 + D_BFx*r_d2*r_dq4*chq2*srq5 + D_BFx*r_d2*r_dq5*chq2*srq5 - D_BFz*r_d2*r_dq4*crq5*shq2 - D_BFz*r_d2*r_dq5*crq5*shq2 + D_BFx*la1*r_dq4*shq2*srq5 + D_BFx*la1*r_dq5*shq2*srq5 + D_BFx*r_d3*r_dq4*shq2*srq5 + D_BFx*r_d3*r_dq5*shq2*srq5 + K_BFx*l3*chq2*crq4*srq5 + K_BFz*l3*chq2*crq5*srq4 - K_BFz*l3*crq4*crq5*shq2 + K_BFx*l3*shq2*srq4*srq5 + K_BFx*la4*chq2*crq5^2*srq4 - K_BFx*la4*crq4*crq5^2*shq2 + K_BFz*la4*chq2*crq5^2*srq4 - K_BFz*la4*crq4*crq5^2*shq2 - D_BFx*l3*r_dq5*chq2*srq4*srq5 + D_BFx*l3*r_dq5*crq4*shq2*srq5 + D_BFz*l3*r_dq5*crq5*shq2*srq4 + K_BFx*la4*chq2*crq4*crq5*srq5 + K_BFz*la4*chq2*crq4*crq5*srq5 + K_BFx*la4*crq5*shq2*srq4*srq5 + K_BFz*la4*crq5*shq2*srq4*srq5 + D_BFz*l3*r_dq5*chq2*crq4*crq5

        ];

%tobeElimitedVars = [];
%expr = eliminate(eqns, tobeElimitedVars);
ddthetaSolved = solve(eqns, [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5]);
%ddthetaSolved = solve(expr, [ddq_2, ddq_3, ddq_4]);


simplify(ddthetaSolved.h_ddq2)
simplify(ddthetaSolved.r_ddd2)
simplify(ddthetaSolved.r_ddd3)
simplify(ddthetaSolved.r_ddq4)
simplify(ddthetaSolved.r_ddq5)