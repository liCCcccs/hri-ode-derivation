
syms l1 l2 l3 l4
syms lc1 la1 lb1 lc2 lc3 lc4 la4 la2
syms q1 h_q2 dq1 h_dq2 ddq1 h_ddq2 
syms r_d2 r_dd2 r_ddd2
syms r_d3 r_dd3 r_ddd3
syms r_q4 r_dq4 r_ddq4
syms r_q5 r_dq5 r_ddq5
syms K_AFz K_AFx K_AMy
syms K_BFz K_BFx K_BMy
syms m1 m2 m3 m4 g
syms I_G1z I_G2z I_G3z I_G4z
syms tau1 tau2 tau3 tau4
q1 = 0; dq1 = 0; ddq1 = 0;

eqns = [tau2 == I_G2z*h_ddq2 + K_BMy*h_q2 - K_BMy*r_q4 - K_BMy*r_q5 + h_ddq2*lc2^2*m2 - K_BFz*la4*lc2*sin(r_q4 - h_q2 + r_q5) + g*lc2*m2*sin(h_q2) + K_BFz*lc2*r_d2*cos(h_q2) - K_BFz*l1*lc2*sin(h_q2) + K_BFz*la1*lc2*sin(h_q2) + K_BFz*lc2*r_d3*sin(h_q2) + K_BFz*l3*lc2*sin(h_q2 - r_q4);

        tau3 == - K_BMy*h_q2 + 2*I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_AMy*r_q4 + K_BMy*r_q4 + K_BMy*r_q5 + l3^2*m4*r_ddq4 + lc3^2*m3*r_ddq4 - K_AFz*lc3*r_d2 + (K_BFx*l3^2*sin(h_q2 - r_q4 + r_q5))/2 + (K_BFx*l3^2*sin(r_q4 - h_q2 + r_q5))/2 - (K_BFz*l3^2*sin(h_q2 - r_q4 + r_q5))/2 + (K_BFz*l3^2*sin(r_q4 - h_q2 + r_q5))/2 + K_BFz*la4*lc4*sin(r_q4 - h_q2 + r_q5) + la4*lc4*m4*r_ddq5 + (K_BFx*l3*r_d2*cos(h_q2 + r_q5))/2 - (K_BFz*l3*r_d2*cos(h_q2 + r_q5))/2 - (K_BFx*l1*l3*sin(h_q2 + r_q5))/2 + (K_BFz*l1*l3*sin(h_q2 + r_q5))/2 + (K_BFx*l3*la1*sin(h_q2 + r_q5))/2 - (K_BFz*l3*la1*sin(h_q2 + r_q5))/2 + (K_BFx*l3*r_d3*sin(h_q2 + r_q5))/2 - (K_BFz*l3*r_d3*sin(h_q2 + r_q5))/2 - lc4*m4*r_ddd2*cos(r_q4 + r_q5) - K_BFz*lc4*r_d2*cos(h_q2) + g*l3*m4*sin(r_q4) + g*lc3*m3*sin(r_q4) + K_BFz*l1*lc4*sin(h_q2) - K_BFz*la1*lc4*sin(h_q2) + (K_BFx*l3*la4*sin(r_q4 - h_q2 + 2*r_q5))/2 + (K_BFz*l3*la4*sin(r_q4 - h_q2 + 2*r_q5))/2 - lc4*m4*r_ddd3*sin(r_q4 + r_q5) - K_BFz*lc4*r_d3*sin(h_q2) - K_BFx*l3*la2*sin(r_q5) - l3*m4*r_ddd2*cos(r_q4) - lc3*m3*r_ddd2*cos(r_q4) - (K_BFx*l3*r_d2*cos(h_q2 - r_q5))/2 - (K_BFz*l3*r_d2*cos(h_q2 - r_q5))/2 + (K_BFx*l1*l3*sin(h_q2 - r_q5))/2 + (K_BFz*l1*l3*sin(h_q2 - r_q5))/2 - (K_BFx*l3*la1*sin(h_q2 - r_q5))/2 + (K_BFx*l3*la4*sin(h_q2 - r_q4))/2 - (K_BFz*l3*la1*sin(h_q2 - r_q5))/2 - (K_BFz*l3*la4*sin(h_q2 - r_q4))/2 - K_BFz*l3*lc4*sin(h_q2 - r_q4) - l3*m4*r_ddd3*sin(r_q4) - lc3*m3*r_ddd3*sin(r_q4) - (K_BFx*l3*r_d3*sin(h_q2 - r_q5))/2 - (K_BFz*l3*r_d3*sin(h_q2 - r_q5))/2 + g*lc4*m4*sin(r_q4 + r_q5) - l3*la4*m4*r_dq5^2*sin(r_q5) + l3*lc4*m4*r_dq4^2*sin(r_q5) + l3*la4*m4*r_ddq5*cos(r_q5) + l3*lc4*m4*r_ddq4*cos(r_q5);
        tau4 == l3*lc4*m4*sin(r_q5)*r_dq4^2 - K_BMy*h_q2 + I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_BMy*r_q4 + K_BMy*r_q5 + K_BFz*la4*lc4*sin(r_q4 - h_q2 + r_q5) + la4*lc4*m4*r_ddq5 - lc4*m4*r_ddd2*cos(r_q4 + r_q5) - K_BFz*lc4*r_d2*cos(h_q2) + K_BFz*l1*lc4*sin(h_q2) - K_BFz*la1*lc4*sin(h_q2) - lc4*m4*r_ddd3*sin(r_q4 + r_q5) - K_BFz*lc4*r_d3*sin(h_q2) - K_BFz*l3*lc4*sin(h_q2 - r_q4) + g*lc4*m4*sin(r_q4 + r_q5) + l3*lc4*m4*r_ddq4*cos(r_q5);

        0 == K_AFx*r_d3 - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 - (K_BFx*l1*cos(h_q2 - r_q5))/2 - (K_BFz*l1*cos(h_q2 - r_q5))/2 + (K_BFx*la1*cos(h_q2 - r_q5))/2 + (K_BFx*la4*cos(h_q2 - r_q4))/2 + (K_BFz*la1*cos(h_q2 - r_q5))/2 - (K_BFz*la4*cos(h_q2 - r_q4))/2 + m3*r_ddd3*cos(r_q4) + m4*r_ddd3*cos(r_q4) + (K_BFx*r_d3*cos(h_q2 - r_q5))/2 + (K_BFz*r_d3*cos(h_q2 - r_q5))/2 - m3*r_ddd2*sin(r_q4) - m4*r_ddd2*sin(r_q4) - (K_BFx*r_d2*sin(h_q2 - r_q5))/2 - (K_BFz*r_d2*sin(h_q2 - r_q5))/2 + (K_BFx*l3*cos(h_q2 - r_q4 + r_q5))/2 + (K_BFx*l3*cos(r_q4 - h_q2 + r_q5))/2 - (K_BFz*l3*cos(h_q2 - r_q4 + r_q5))/2 + (K_BFz*l3*cos(r_q4 - h_q2 + r_q5))/2 - (K_BFx*l1*cos(h_q2 + r_q5))/2 + (K_BFz*l1*cos(h_q2 + r_q5))/2 + (K_BFx*la1*cos(h_q2 + r_q5))/2 - (K_BFz*la1*cos(h_q2 + r_q5))/2 + (K_BFx*r_d3*cos(h_q2 + r_q5))/2 - (K_BFz*r_d3*cos(h_q2 + r_q5))/2 - (K_BFx*r_d2*sin(h_q2 + r_q5))/2 + (K_BFz*r_d2*sin(h_q2 + r_q5))/2 - g*m3*cos(r_q4) - g*m4*cos(r_q4) + (K_BFx*la4*cos(r_q4 - h_q2 + 2*r_q5))/2 + (K_BFz*la4*cos(r_q4 - h_q2 + 2*r_q5))/2 - K_BFx*la2*cos(r_q5) - la4*m4*r_ddq5*sin(r_q5) - la4*m4*r_dq5^2*cos(r_q5);
        0 == (K_BFx*l1*sin(h_q2 - r_q5))/2 - m3*r_ddd2*cos(r_q4) - m4*r_ddd2*cos(r_q4) - (K_BFx*r_d2*cos(h_q2 - r_q5))/2 - (K_BFz*r_d2*cos(h_q2 - r_q5))/2 - K_AFz*r_d2 + (K_BFz*l1*sin(h_q2 - r_q5))/2 - (K_BFx*la1*sin(h_q2 - r_q5))/2 + (K_BFx*la4*sin(h_q2 - r_q4))/2 - (K_BFz*la1*sin(h_q2 - r_q5))/2 - (K_BFz*la4*sin(h_q2 - r_q4))/2 - m3*r_ddd3*sin(r_q4) - m4*r_ddd3*sin(r_q4) - (K_BFx*r_d3*sin(h_q2 - r_q5))/2 - (K_BFz*r_d3*sin(h_q2 - r_q5))/2 + (K_BFx*l3*sin(h_q2 - r_q4 + r_q5))/2 + (K_BFx*l3*sin(r_q4 - h_q2 + r_q5))/2 - (K_BFz*l3*sin(h_q2 - r_q4 + r_q5))/2 + (K_BFz*l3*sin(r_q4 - h_q2 + r_q5))/2 + l3*m4*r_ddq4 + lc3*m3*r_ddq4 + (K_BFx*r_d2*cos(h_q2 + r_q5))/2 - (K_BFz*r_d2*cos(h_q2 + r_q5))/2 - (K_BFx*l1*sin(h_q2 + r_q5))/2 + (K_BFz*l1*sin(h_q2 + r_q5))/2 + (K_BFx*la1*sin(h_q2 + r_q5))/2 - (K_BFz*la1*sin(h_q2 + r_q5))/2 + (K_BFx*r_d3*sin(h_q2 + r_q5))/2 - (K_BFz*r_d3*sin(h_q2 + r_q5))/2 + g*m3*sin(r_q4) + g*m4*sin(r_q4) + (K_BFx*la4*sin(r_q4 - h_q2 + 2*r_q5))/2 + (K_BFz*la4*sin(r_q4 - h_q2 + 2*r_q5))/2 - K_BFx*la2*sin(r_q5) + la4*m4*r_ddq5*cos(r_q5) - la4*m4*r_dq5^2*sin(r_q5)

        ];

ddthetaSolved = solve(eqns, [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5]);
