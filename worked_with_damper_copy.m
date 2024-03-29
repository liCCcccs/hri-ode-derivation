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
q1 = pi/2; dq1 = 0; ddq1 = 0;
% D_AFz = 0; D_AFx = 0; D_AMy = 0;
% D_BFz = 0; D_BFx = 0; D_BMy = 0;
%K_AMy=0;
%K_BMy=0;
tau1 = 0; tau3 = 0;
%h_q2 = 0; h_dq2 = 0; h_ddq2 = 0;
%r_q4 = 0; r_dq4 = 0; r_ddq4 = 0;
%r_q5 = 0; r_dq5 = 0; r_ddq5 = 0;
%r_d2 = 0; r_dd2 = 0; r_ddd2 = 0;
%r_d3 = 0; r_dd3 = 0; r_ddd3 = 0; 

% intermidiate variables
syms f_intBi  % D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))))
syms f_intBj  % - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))))
syms m4gi  % m4*g*cos(q1+r_q4+r_q5)
syms m4gj  % -m4*g*sin(q1+r_q4+r_q5)
syms tau_intBz  % K_BMy*(r_q4 + r_q5 - h_q2) + D_BMy*(r_dq4 + r_dq5 - h_dq2)
syms f_intAi  % K_AFx * r_d3 + D_AFx * r_dd3
syms f_intAj  % + K_AFz * r_d2 + D_AFz * r_dd2
syms m3gi  % m3*g*cos(q1+r_q4)
syms m3gj  % -m3*g*sin(q1+r_q4)
syms tau_intAz  % K_AMy*r_q4 + D_AMy*r_dq4
syms m2gi  % m2*g*cos(q1+h_q2)
syms m2gj  % -m2*g*sin(q1+h_q2)
syms m1gi  % m1*g*cos(q1)
syms m1gj  % -m1*g*sin(q1)

% eqns = [%tau1 = D_BMy*h_dq2 + I_G1z*ddq1 + I_G2z*ddq1 + I_G2z*h_ddq2 + K_BMy*h_q2 + D_AMy*r_dq4 - D_BMy*r_dq4 - D_BMy*r_dq5 + K_AMy*r_q4 - K_BMy*r_q4 - K_BMy*r_q5 + h_ddq2*lc2^2*m2 - (D_BFx*l1*r_dd2)/2 - (D_BFx*l1*r_dd3)/2 - (D_BFz*l1*r_dd2)/2 - (D_BFz*l1*r_dd3)/2 + D_AFz*lc1*r_dd2 + D_BFx*lc1*r_dd2 + D_BFx*lc1*r_dd3 + D_BFz*lc1*r_dd2 + D_BFz*lc1*r_dd3 - (K_BFx*l1*r_d2)/2 - (K_BFz*l1*r_d2)/2 + K_AFz*lc1*r_d2 + K_BFx*lc1*r_d2 + K_BFz*lc1*r_d2 - (K_BFx*l1^2*sin(2*h_q2))/2 + (K_BFz*l1^2*sin(2*h_q2))/2 - (D_BFx*h_dq2*l1^2)/2 - (D_BFz*h_dq2*l1^2)/2 + ddq1*lc1^2*m1 + ddq1*lc2^2*m2 - (K_BFx*l1*la4*sin(r_q4 - 2*h_q2 + r_q5))/2 + (K_BFz*l1*la4*sin(r_q4 - 2*h_q2 + r_q5))/2 + (D_BFx*h_dq2*l1^2*cos(2*h_q2))/2 - (D_BFz*h_dq2*l1^2*cos(2*h_q2))/2 + K_BFx*la4*lc1*sin(r_q4 - 2*h_q2 + r_q5) - K_BFz*la4*lc1*sin(r_q4 - 2*h_q2 + r_q5) - K_BFz*la4*lc2*sin(r_q4 - h_q2 + r_q5) + (ddq1*l1*lc1*m2)/2 + (h_ddq2*l1*lc1*m2)/2 - l1*m4*r_ddd2*cos(r_q4 - h_q2 + r_q5) + lc1*m4*r_ddd2*cos(r_q4 - h_q2 + r_q5) - l1*m4*r_ddd3*sin(r_q4 - h_q2 + r_q5) + lc1*m4*r_ddd3*sin(r_q4 - h_q2 + r_q5) + (K_BFx*l1*la4*sin(r_q4 + r_q5))/2 + (K_BFz*l1*la4*sin(r_q4 + r_q5))/2 - K_BFx*la4*lc1*sin(r_q4 + r_q5) - K_BFz*la4*lc1*sin(r_q4 + r_q5) + g*lc2*m2*sin(h_q2 + q1) + D_BFz*lc2*r_dd2*cos(h_q2) + D_BFz*lc2*r_dd3*cos(h_q2) + K_BFz*lc2*r_d2*cos(h_q2) - K_BFx*l1*la2*sin(h_q2) - K_BFz*l1*lc2*sin(h_q2) + 2*K_BFx*la2*lc1*sin(h_q2) + K_BFz*la1*lc2*sin(h_q2) + K_BFz*lc2*r_d3*sin(h_q2) + (K_BFx*l1*l3*sin(r_q4))/2 + (K_BFz*l1*l3*sin(r_q4))/2 - K_BFx*l3*lc1*sin(r_q4) - K_BFz*l3*lc1*sin(r_q4) + g*lc1*m1*(-1) + g*lc1*m2*(-1) + g*l1*m4*sin(q1 - h_q2 + r_q4 + r_q5) + K_BFz*l3*lc2*sin(h_q2 - r_q4) - g*lc1*m4*sin(q1 - h_q2 + r_q4 + r_q5) + (D_BFx*l1*r_dd2*cos(2*h_q2))/2 + (D_BFx*l1*r_dd3*cos(2*h_q2))/2 - (D_BFz*l1*r_dd2*cos(2*h_q2))/2 - (D_BFz*l1*r_dd3*cos(2*h_q2))/2 - D_BFx*lc1*r_dd2*cos(2*h_q2) - D_BFx*lc1*r_dd3*cos(2*h_q2) + D_BFz*lc1*r_dd2*cos(2*h_q2) + D_BFz*lc1*r_dd3*cos(2*h_q2) + (K_BFx*l1*r_d2*cos(2*h_q2))/2 - (K_BFz*l1*r_d2*cos(2*h_q2))/2 - K_BFx*lc1*r_d2*cos(2*h_q2) + K_BFz*lc1*r_d2*cos(2*h_q2) + (K_BFx*l1*la1*sin(2*h_q2))/2 - (K_BFz*l1*la1*sin(2*h_q2))/2 + K_BFx*l1*lc1*sin(2*h_q2) - K_BFz*l1*lc1*sin(2*h_q2) - K_BFx*la1*lc1*sin(2*h_q2) + K_BFz*la1*lc1*sin(2*h_q2) + (K_BFx*l1*r_d3*sin(2*h_q2))/2 - (K_BFz*l1*r_d3*sin(2*h_q2))/2 - K_BFx*lc1*r_d3*sin(2*h_q2) + K_BFz*lc1*r_d3*sin(2*h_q2) + D_BFx*h_dq2*l1*lc1 + D_BFz*h_dq2*l1*lc1 + (D_BFx*l1*la1*r_dq4)/2 + (D_BFx*l1*la1*r_dq5)/2 + (D_BFz*l1*la1*r_dq4)/2 + (D_BFz*l1*la1*r_dq5)/2 - D_BFx*la1*lc1*r_dq4 - D_BFx*la1*lc1*r_dq5 - D_BFz*la1*lc1*r_dq4 - D_BFz*la1*lc1*r_dq5 + (K_BFx*l1*l3*sin(2*h_q2 - r_q4))/2 - (K_BFz*l1*l3*sin(2*h_q2 - r_q4))/2 - K_BFx*l3*lc1*sin(2*h_q2 - r_q4) + K_BFz*l3*lc1*sin(2*h_q2 - r_q4) + (D_BFx*l1*r_d3*r_dq4)/2 + (D_BFx*l1*r_d3*r_dq5)/2 + (D_BFz*l1*r_d3*r_dq4)/2 + (D_BFz*l1*r_d3*r_dq5)/2 - D_BFx*lc1*r_d3*r_dq4 - D_BFx*lc1*r_d3*r_dq5 - D_BFz*lc1*r_d3*r_dq4 - D_BFz*lc1*r_d3*r_dq5 - D_BFx*h_dq2*l1*lc1*cos(2*h_q2) + D_BFz*h_dq2*l1*lc1*cos(2*h_q2) + l1*l3*m4*r_ddq4*cos(h_q2 - r_q5) - l3*lc1*m4*r_ddq4*cos(h_q2 - r_q5) - (D_BFx*l1*la1*r_dq4*cos(2*h_q2))/2 - (D_BFx*l1*la1*r_dq5*cos(2*h_q2))/2 + (D_BFz*l1*la1*r_dq4*cos(2*h_q2))/2 + (D_BFz*l1*la1*r_dq5*cos(2*h_q2))/2 + D_BFx*la1*lc1*r_dq4*cos(2*h_q2) + D_BFx*la1*lc1*r_dq5*cos(2*h_q2) - D_BFz*la1*lc1*r_dq4*cos(2*h_q2) - D_BFz*la1*lc1*r_dq5*cos(2*h_q2) - (D_BFx*l1*r_d3*r_dq4*cos(2*h_q2))/2 - (D_BFx*l1*r_d3*r_dq5*cos(2*h_q2))/2 + (D_BFz*l1*r_d3*r_dq4*cos(2*h_q2))/2 + (D_BFz*l1*r_d3*r_dq5*cos(2*h_q2))/2 + D_BFx*lc1*r_d3*r_dq4*cos(2*h_q2) + D_BFx*lc1*r_d3*r_dq5*cos(2*h_q2) - D_BFz*lc1*r_d3*r_dq4*cos(2*h_q2) - D_BFz*lc1*r_d3*r_dq5*cos(2*h_q2) + (D_BFx*l1*r_d2*r_dq4*sin(2*h_q2))/2 + (D_BFx*l1*r_d2*r_dq5*sin(2*h_q2))/2 - (D_BFz*l1*r_d2*r_dq4*sin(2*h_q2))/2 - (D_BFz*l1*r_d2*r_dq5*sin(2*h_q2))/2 - D_BFx*lc1*r_d2*r_dq4*sin(2*h_q2) - D_BFx*lc1*r_d2*r_dq5*sin(2*h_q2) + D_BFz*lc1*r_d2*r_dq4*sin(2*h_q2) + D_BFz*lc1*r_d2*r_dq5*sin(2*h_q2) + (ddq1*l1*lc1*m2*cos(2*h_q2))/2 - (h_ddq2*l1*lc1*m2*cos(2*h_q2))/2 + dq1^2*l1*lc2*m2*sin(h_q2) - dq1^2*lc1*lc2*m2*sin(h_q2) - h_dq2^2*lc1*lc2*m2*sin(h_q2) - (D_BFx*l1*l3*r_dq5*cos(2*h_q2 - r_q4))/2 + (D_BFz*l1*l3*r_dq5*cos(2*h_q2 - r_q4))/2 + D_BFx*l3*lc1*r_dq5*cos(2*h_q2 - r_q4) - D_BFz*l3*lc1*r_dq5*cos(2*h_q2 - r_q4) - l1*la4*m4*r_dq5^2*sin(h_q2) + la4*lc1*m4*r_dq5^2*sin(h_q2) - l1*l3*m4*r_dq4^2*sin(h_q2 - r_q5) + l3*lc1*m4*r_dq4^2*sin(h_q2 - r_q5) + ddq1*l1*la1*m4*sin(r_q4 - h_q2 + r_q5) - ddq1*la1*lc1*m4*sin(r_q4 - h_q2 + r_q5) + D_BFz*h_dq2*l1*lc2*cos(h_q2) - D_BFz*la1*lc2*r_dq4*cos(h_q2) - D_BFz*la1*lc2*r_dq5*cos(h_q2) - D_BFz*lc2*r_d3*r_dq4*cos(h_q2) - D_BFz*lc2*r_d3*r_dq5*cos(h_q2) + (D_BFx*l1*l3*r_dq5*cos(r_q4))/2 + (D_BFz*l1*l3*r_dq5*cos(r_q4))/2 - D_BFx*l3*lc1*r_dq5*cos(r_q4) - D_BFz*l3*lc1*r_dq5*cos(r_q4) + D_BFz*lc2*r_d2*r_dq4*sin(h_q2) + D_BFz*lc2*r_d2*r_dq5*sin(h_q2) + ddq1*l1*lc2*m2*cos(h_q2) + ddq1*lc1*lc2*m2*cos(h_q2) + dq1^2*l1*la1*m4*cos(r_q4 - h_q2 + r_q5) - dq1^2*la1*lc1*m4*cos(r_q4 - h_q2 + r_q5) + h_ddq2*lc1*lc2*m2*cos(h_q2) + l1*la4*m4*r_ddq5*cos(h_q2) - la4*lc1*m4*r_ddq5*cos(h_q2) - D_BFz*l3*lc2*r_dq5*cos(h_q2 - r_q4) - 2*dq1*h_dq2*lc1*lc2*m2*sin(h_q2)
% 
%         tau2 == D_BMy*h_dq2 + I_G2z*ddq1 + I_G2z*h_ddq2 + K_BMy*h_q2 - D_BMy*r_dq4 - D_BMy*r_dq5 - K_BMy*r_q4 - K_BMy*r_q5 + h_ddq2*lc2^2*m2 + ddq1*lc2^2*m2 - K_BFz*la4*lc2*sin(r_q4 - h_q2 + r_q5) + g*lc2*m2*sin(h_q2 + q1) + K_BFz*lc2*r_d3*cos(h_q2) + D_BFz*lc2*r_dd2*sin(h_q2) + D_BFz*lc2*r_dd3*sin(h_q2) - K_BFz*l1*lc2*sin(h_q2) + K_BFz*la1*lc2*sin(h_q2) + K_BFz*lc2*r_d2*sin(h_q2) + K_BFz*l3*lc2*sin(h_q2 - r_q4) + dq1^2*l1*lc2*m2*sin(h_q2) + D_BFz*h_dq2*l1*lc2*cos(h_q2) - D_BFz*la1*lc2*r_dq4*cos(h_q2) - D_BFz*la1*lc2*r_dq5*cos(h_q2) - D_BFz*lc2*r_d2*r_dq4*cos(h_q2) - D_BFz*lc2*r_d2*r_dq5*cos(h_q2) + D_BFz*lc2*r_d3*r_dq4*sin(h_q2) + D_BFz*lc2*r_d3*r_dq5*sin(h_q2) + ddq1*l1*lc2*m2*cos(h_q2) - D_BFz*l3*lc2*r_dq5*cos(h_q2 - r_q4)
% 
%         tau3 == 2*I_G4z*ddq1 - D_BMy*h_dq2 - K_BMy*h_q2 + D_AMy*r_dq4 + D_BMy*r_dq4 + D_BMy*r_dq5 + 2*I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_AMy*r_q4 + K_BMy*r_q4 + K_BMy*r_q5 + l3^2*m4*r_ddq4 + lc3^2*m3*r_ddq4 - D_AFz*lc3*r_dd2 - K_AFz*lc3*r_d2 + la4*lc4*m4*r_ddq5 - K_BFz*lc4*r_d3*cos(h_q2) - D_BFz*lc4*r_dd2*sin(h_q2) - D_BFz*lc4*r_dd3*sin(h_q2) + K_BFz*l1*lc4*sin(h_q2) - K_BFz*la1*lc4*sin(h_q2) - K_BFz*lc4*r_d2*sin(h_q2) - K_BFx*l3*la2*sin(r_q5) - l3*m4*r_ddd3*cos(r_q4) - lc3*m3*r_ddd3*cos(r_q4) - l3*m4*r_ddd2*sin(r_q4) - lc3*m3*r_ddd2*sin(r_q4) + K_BFx*l3^2*sin(h_q2)*sin(r_q4)*sin(r_q5) - K_BFz*l3*r_d3*cos(h_q2)*cos(r_q5) + D_BFx*l3*r_dd2*cos(h_q2)*sin(r_q5) + D_BFx*l3*r_dd3*cos(h_q2)*sin(r_q5) - D_BFz*l3*r_dd2*cos(r_q5)*sin(h_q2) - D_BFz*l3*r_dd3*cos(r_q5)*sin(h_q2) - K_BFx*l1*l3*cos(h_q2)*sin(r_q5) + K_BFz*l1*l3*cos(r_q5)*sin(h_q2) + K_BFx*l3*la1*cos(h_q2)*sin(r_q5) - K_BFx*l3*la4*cos(h_q2)*sin(r_q4) + K_BFx*l3*la4*cos(r_q4)*sin(h_q2) - K_BFz*l3*la1*cos(r_q5)*sin(h_q2) + K_BFz*l3*lc4*cos(h_q2)*sin(r_q4) - K_BFz*l3*lc4*cos(r_q4)*sin(h_q2) + K_BFx*l3*r_d2*cos(h_q2)*sin(r_q5) - K_BFz*l3*r_d2*cos(r_q5)*sin(h_q2) - K_BFx*l3*r_d3*sin(h_q2)*sin(r_q5) + dq1^2*l3*la1*m4*sin(r_q4) + dq1^2*la1*lc3*m3*sin(r_q4) - lc4*m4*r_ddd3*cos(r_q4)*cos(r_q5) + g*l3*m4*cos(q1)*sin(r_q4) + g*l3*m4*cos(r_q4)*sin(q1) + g*lc3*m3*cos(q1)*sin(r_q4) + g*lc3*m3*cos(r_q4)*sin(q1) - l3*la4*m4*r_dq5^2*sin(r_q5) + l3*lc4*m4*r_dq4^2*sin(r_q5) - lc4*m4*r_ddd2*cos(r_q4)*sin(r_q5) - lc4*m4*r_ddd2*cos(r_q5)*sin(r_q4) + lc4*m4*r_ddd3*sin(r_q4)*sin(r_q5) - D_BFz*h_dq2*l1*lc4*cos(h_q2) + D_BFz*la1*lc4*r_dq4*cos(h_q2) + D_BFz*la1*lc4*r_dq5*cos(h_q2) + D_BFz*lc4*r_d2*r_dq4*cos(h_q2) + D_BFz*lc4*r_d2*r_dq5*cos(h_q2) - D_BFz*lc4*r_d3*r_dq4*sin(h_q2) - D_BFz*lc4*r_d3*r_dq5*sin(h_q2) - ddq1*l3*la1*m4*cos(r_q4) - ddq1*la1*lc3*m3*cos(r_q4) + l3*la4*m4*r_ddq5*cos(r_q5) + l3*lc4*m4*r_ddq4*cos(r_q5) + K_BFx*l3^2*cos(h_q2)*cos(r_q4)*sin(r_q5) + K_BFz*l3^2*cos(h_q2)*cos(r_q5)*sin(r_q4) - K_BFz*l3^2*cos(r_q4)*cos(r_q5)*sin(h_q2) + K_BFz*la4*lc4*cos(h_q2)*cos(r_q4)*sin(r_q5) + K_BFz*la4*lc4*cos(h_q2)*cos(r_q5)*sin(r_q4) - K_BFz*la4*lc4*cos(r_q4)*cos(r_q5)*sin(h_q2) + dq1^2*la1*lc4*m4*cos(r_q4)*sin(r_q5) + dq1^2*la1*lc4*m4*cos(r_q5)*sin(r_q4) + K_BFz*la4*lc4*sin(h_q2)*sin(r_q4)*sin(r_q5) + g*lc4*m4*cos(q1)*cos(r_q4)*sin(r_q5) + g*lc4*m4*cos(q1)*cos(r_q5)*sin(r_q4) + g*lc4*m4*cos(r_q4)*cos(r_q5)*sin(q1) - g*lc4*m4*sin(q1)*sin(r_q4)*sin(r_q5) + D_BFz*l3^2*r_dq5*cos(h_q2)*cos(r_q4)*cos(r_q5) + K_BFx*l3*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFx*l3*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) + K_BFz*l3*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFz*l3*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) - D_BFx*l3^2*r_dq5*cos(h_q2)*sin(r_q4)*sin(r_q5) + D_BFx*l3^2*r_dq5*cos(r_q4)*sin(h_q2)*sin(r_q5) + D_BFz*l3^2*r_dq5*cos(r_q5)*sin(h_q2)*sin(r_q4) - D_BFz*h_dq2*l1*l3*cos(h_q2)*cos(r_q5) + D_BFz*l3*la1*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*l3*la1*r_dq5*cos(h_q2)*cos(r_q5) + D_BFz*l3*lc4*r_dq5*cos(h_q2)*cos(r_q4) + D_BFz*l3*r_d2*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*l3*r_d2*r_dq5*cos(h_q2)*cos(r_q5) - D_BFx*h_dq2*l1*l3*sin(h_q2)*sin(r_q5) + D_BFx*l3*r_d3*r_dq4*cos(h_q2)*sin(r_q5) + D_BFx*l3*r_d3*r_dq5*cos(h_q2)*sin(r_q5) - D_BFz*l3*r_d3*r_dq4*cos(r_q5)*sin(h_q2) - D_BFz*l3*r_d3*r_dq5*cos(r_q5)*sin(h_q2) + D_BFx*l3*la1*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*l3*la1*r_dq5*sin(h_q2)*sin(r_q5) + D_BFz*l3*lc4*r_dq5*sin(h_q2)*sin(r_q4) - ddq1*la1*lc4*m4*cos(r_q4)*cos(r_q5) + D_BFx*l3*r_d2*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*l3*r_d2*r_dq5*sin(h_q2)*sin(r_q5) + ddq1*la1*lc4*m4*sin(r_q4)*sin(r_q5) + K_BFx*l3*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFz*l3*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFx*l3*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5) + K_BFz*l3*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5)
% 
%         tau4 == I_G4z*ddq1 - D_BMy*h_dq2 - K_BMy*h_q2 + D_BMy*r_dq4 + D_BMy*r_dq5 + I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_BMy*r_q4 + K_BMy*r_q5 + K_BFz*la4*lc4*sin(r_q4 - h_q2 + r_q5) + la4*lc4*m4*r_ddq5 - lc4*m4*r_ddd3*cos(r_q4 + r_q5) - K_BFz*lc4*r_d3*cos(h_q2) - D_BFz*lc4*r_dd2*sin(h_q2) - D_BFz*lc4*r_dd3*sin(h_q2) + K_BFz*l1*lc4*sin(h_q2) - K_BFz*la1*lc4*sin(h_q2) - lc4*m4*r_ddd2*sin(r_q4 + r_q5) - K_BFz*lc4*r_d2*sin(h_q2) - K_BFz*l3*lc4*sin(h_q2 - r_q4) + g*lc4*m4*sin(q1 + r_q4 + r_q5) + dq1^2*la1*lc4*m4*sin(r_q4 + r_q5) + l3*lc4*m4*r_dq4^2*sin(r_q5) - D_BFz*h_dq2*l1*lc4*cos(h_q2) - ddq1*la1*lc4*m4*cos(r_q4 + r_q5) + D_BFz*la1*lc4*r_dq4*cos(h_q2) + D_BFz*la1*lc4*r_dq5*cos(h_q2) + D_BFz*lc4*r_d2*r_dq4*cos(h_q2) + D_BFz*lc4*r_d2*r_dq5*cos(h_q2) - D_BFz*lc4*r_d3*r_dq4*sin(h_q2) - D_BFz*lc4*r_d3*r_dq5*sin(h_q2) + D_BFz*l3*lc4*r_dq5*cos(h_q2 - r_q4) + l3*lc4*m4*r_ddq4*cos(r_q5)
% 
%         0 == D_AFx*r_dd3 + K_AFx*r_d3 - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 + m3*r_ddd2*cos(r_q4) + m4*r_ddd2*cos(r_q4) - m3*r_ddd3*sin(r_q4) - m4*r_ddd3*sin(r_q4) - K_BFx*la2*cos(r_q5) - ddq1*la1*m3*sin(r_q4) - ddq1*la1*m4*sin(r_q4) - la4*m4*r_ddq5*sin(r_q5) + D_BFx*r_dd2*cos(h_q2)*cos(r_q5) + D_BFx*r_dd3*cos(h_q2)*cos(r_q5) - K_BFx*l1*cos(h_q2)*cos(r_q5) + K_BFx*la1*cos(h_q2)*cos(r_q5) - K_BFz*la4*cos(h_q2)*cos(r_q4) + K_BFx*r_d2*cos(h_q2)*cos(r_q5) - K_BFx*r_d3*cos(r_q5)*sin(h_q2) + K_BFz*r_d3*cos(h_q2)*sin(r_q5) - dq1^2*la1*m3*cos(r_q4) - dq1^2*la1*m4*cos(r_q4) + D_BFz*r_dd2*sin(h_q2)*sin(r_q5) + D_BFz*r_dd3*sin(h_q2)*sin(r_q5) - K_BFz*l1*sin(h_q2)*sin(r_q5) + K_BFz*la1*sin(h_q2)*sin(r_q5) - K_BFz*la4*sin(h_q2)*sin(r_q4) - g*m3*cos(q1)*cos(r_q4) - g*m4*cos(q1)*cos(r_q4) + K_BFz*r_d2*sin(h_q2)*sin(r_q5) - la4*m4*r_dq5^2*cos(r_q5) + g*m3*sin(q1)*sin(r_q4) + g*m4*sin(q1)*sin(r_q4) + K_BFx*la4*cos(r_q5)^2*sin(h_q2)*sin(r_q4) + K_BFz*la4*cos(r_q5)^2*sin(h_q2)*sin(r_q4) - D_BFx*h_dq2*l1*cos(r_q5)*sin(h_q2) + D_BFz*h_dq2*l1*cos(h_q2)*sin(r_q5) + D_BFx*r_d3*r_dq4*cos(h_q2)*cos(r_q5) + D_BFx*r_d3*r_dq5*cos(h_q2)*cos(r_q5) + D_BFx*la1*r_dq4*cos(r_q5)*sin(h_q2) + D_BFx*la1*r_dq5*cos(r_q5)*sin(h_q2) - D_BFz*la1*r_dq4*cos(h_q2)*sin(r_q5) - D_BFz*la1*r_dq5*cos(h_q2)*sin(r_q5) + D_BFx*r_d2*r_dq4*cos(r_q5)*sin(h_q2) + D_BFx*r_d2*r_dq5*cos(r_q5)*sin(h_q2) - D_BFz*r_d2*r_dq4*cos(h_q2)*sin(r_q5) - D_BFz*r_d2*r_dq5*cos(h_q2)*sin(r_q5) + D_BFz*r_d3*r_dq4*sin(h_q2)*sin(r_q5) + D_BFz*r_d3*r_dq5*sin(h_q2)*sin(r_q5) + K_BFx*l3*cos(h_q2)*cos(r_q4)*cos(r_q5) + K_BFx*l3*cos(r_q5)*sin(h_q2)*sin(r_q4) - K_BFz*l3*cos(h_q2)*sin(r_q4)*sin(r_q5) + K_BFz*l3*cos(r_q4)*sin(h_q2)*sin(r_q5) + K_BFx*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)^2 + K_BFz*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)^2 - D_BFx*l3*r_dq5*cos(h_q2)*cos(r_q5)*sin(r_q4) + D_BFx*l3*r_dq5*cos(r_q4)*cos(r_q5)*sin(h_q2) - D_BFz*l3*r_dq5*cos(h_q2)*cos(r_q4)*sin(r_q5) - D_BFz*l3*r_dq5*sin(h_q2)*sin(r_q4)*sin(r_q5) - K_BFx*la4*cos(h_q2)*cos(r_q5)*sin(r_q4)*sin(r_q5) + K_BFx*la4*cos(r_q4)*cos(r_q5)*sin(h_q2)*sin(r_q5) - K_BFz*la4*cos(h_q2)*cos(r_q5)*sin(r_q4)*sin(r_q5) + K_BFz*la4*cos(r_q4)*cos(r_q5)*sin(h_q2)*sin(r_q5)
%         0 == l3*m4*r_ddq4 - K_AFz*r_d2 - m3*r_ddd3*cos(r_q4) - m4*r_ddd3*cos(r_q4) - m3*r_ddd2*sin(r_q4) - m4*r_ddd2*sin(r_q4) - D_AFz*r_dd2 + lc3*m3*r_ddq4 - K_BFx*la2*sin(r_q5) - ddq1*la1*m3*cos(r_q4) - ddq1*la1*m4*cos(r_q4) + la4*m4*r_ddq5*cos(r_q5) - K_BFz*r_d3*cos(h_q2)*cos(r_q5) + D_BFx*r_dd2*cos(h_q2)*sin(r_q5) + D_BFx*r_dd3*cos(h_q2)*sin(r_q5) - D_BFz*r_dd2*cos(r_q5)*sin(h_q2) - D_BFz*r_dd3*cos(r_q5)*sin(h_q2) - K_BFx*l1*cos(h_q2)*sin(r_q5) + K_BFz*l1*cos(r_q5)*sin(h_q2) + K_BFx*la1*cos(h_q2)*sin(r_q5) - K_BFx*la4*cos(h_q2)*sin(r_q4) + K_BFx*la4*cos(r_q4)*sin(h_q2) - K_BFz*la1*cos(r_q5)*sin(h_q2) + K_BFx*r_d2*cos(h_q2)*sin(r_q5) - K_BFz*r_d2*cos(r_q5)*sin(h_q2) - K_BFx*r_d3*sin(h_q2)*sin(r_q5) + dq1^2*la1*m3*sin(r_q4) + dq1^2*la1*m4*sin(r_q4) + g*m3*cos(q1)*sin(r_q4) + g*m3*cos(r_q4)*sin(q1) + g*m4*cos(q1)*sin(r_q4) + g*m4*cos(r_q4)*sin(q1) - la4*m4*r_dq5^2*sin(r_q5) - D_BFz*h_dq2*l1*cos(h_q2)*cos(r_q5) + D_BFz*la1*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*la1*r_dq5*cos(h_q2)*cos(r_q5) + D_BFz*r_d2*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*r_d2*r_dq5*cos(h_q2)*cos(r_q5) - D_BFx*h_dq2*l1*sin(h_q2)*sin(r_q5) + D_BFx*r_d3*r_dq4*cos(h_q2)*sin(r_q5) + D_BFx*r_d3*r_dq5*cos(h_q2)*sin(r_q5) - D_BFz*r_d3*r_dq4*cos(r_q5)*sin(h_q2) - D_BFz*r_d3*r_dq5*cos(r_q5)*sin(h_q2) + D_BFx*la1*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*la1*r_dq5*sin(h_q2)*sin(r_q5) + D_BFx*r_d2*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*r_d2*r_dq5*sin(h_q2)*sin(r_q5) + K_BFx*l3*cos(h_q2)*cos(r_q4)*sin(r_q5) + K_BFz*l3*cos(h_q2)*cos(r_q5)*sin(r_q4) - K_BFz*l3*cos(r_q4)*cos(r_q5)*sin(h_q2) + K_BFx*l3*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFx*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFx*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) + K_BFz*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFz*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) - D_BFx*l3*r_dq5*cos(h_q2)*sin(r_q4)*sin(r_q5) + D_BFx*l3*r_dq5*cos(r_q4)*sin(h_q2)*sin(r_q5) + D_BFz*l3*r_dq5*cos(r_q5)*sin(h_q2)*sin(r_q4) + K_BFx*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5) + K_BFz*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5) + K_BFx*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFz*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + D_BFz*l3*r_dq5*cos(h_q2)*cos(r_q4)*cos(r_q5)
% 
%         ];

% Manually assigned q1 = pi/2 to avoid numerical error
% eqns = [%tau1 = D_BMy*h_dq2 + I_G1z*ddq1 + I_G2z*ddq1 + I_G2z*h_ddq2 + K_BMy*h_q2 + D_AMy*r_dq4 - D_BMy*r_dq4 - D_BMy*r_dq5 + K_AMy*r_q4 - K_BMy*r_q4 - K_BMy*r_q5 + h_ddq2*lc2^2*m2 - (D_BFx*l1*r_dd2)/2 - (D_BFx*l1*r_dd3)/2 - (D_BFz*l1*r_dd2)/2 - (D_BFz*l1*r_dd3)/2 + D_AFz*lc1*r_dd2 + D_BFx*lc1*r_dd2 + D_BFx*lc1*r_dd3 + D_BFz*lc1*r_dd2 + D_BFz*lc1*r_dd3 - (K_BFx*l1*r_d2)/2 - (K_BFz*l1*r_d2)/2 + K_AFz*lc1*r_d2 + K_BFx*lc1*r_d2 + K_BFz*lc1*r_d2 - (K_BFx*l1^2*sin(2*h_q2))/2 + (K_BFz*l1^2*sin(2*h_q2))/2 - (D_BFx*h_dq2*l1^2)/2 - (D_BFz*h_dq2*l1^2)/2 + ddq1*lc1^2*m1 + ddq1*lc2^2*m2 - (K_BFx*l1*la4*sin(r_q4 - 2*h_q2 + r_q5))/2 + (K_BFz*l1*la4*sin(r_q4 - 2*h_q2 + r_q5))/2 + (D_BFx*h_dq2*l1^2*cos(2*h_q2))/2 - (D_BFz*h_dq2*l1^2*cos(2*h_q2))/2 + K_BFx*la4*lc1*sin(r_q4 - 2*h_q2 + r_q5) - K_BFz*la4*lc1*sin(r_q4 - 2*h_q2 + r_q5) - K_BFz*la4*lc2*sin(r_q4 - h_q2 + r_q5) + (ddq1*l1*lc1*m2)/2 + (h_ddq2*l1*lc1*m2)/2 - l1*m4*r_ddd2*cos(r_q4 - h_q2 + r_q5) + lc1*m4*r_ddd2*cos(r_q4 - h_q2 + r_q5) - l1*m4*r_ddd3*sin(r_q4 - h_q2 + r_q5) + lc1*m4*r_ddd3*sin(r_q4 - h_q2 + r_q5) + (K_BFx*l1*la4*sin(r_q4 + r_q5))/2 + (K_BFz*l1*la4*sin(r_q4 + r_q5))/2 - K_BFx*la4*lc1*sin(r_q4 + r_q5) - K_BFz*la4*lc1*sin(r_q4 + r_q5) + g*lc2*m2*sin(h_q2 + q1) + D_BFz*lc2*r_dd2*cos(h_q2) + D_BFz*lc2*r_dd3*cos(h_q2) + K_BFz*lc2*r_d2*cos(h_q2) - K_BFx*l1*la2*sin(h_q2) - K_BFz*l1*lc2*sin(h_q2) + 2*K_BFx*la2*lc1*sin(h_q2) + K_BFz*la1*lc2*sin(h_q2) + K_BFz*lc2*r_d3*sin(h_q2) + (K_BFx*l1*l3*sin(r_q4))/2 + (K_BFz*l1*l3*sin(r_q4))/2 - K_BFx*l3*lc1*sin(r_q4) - K_BFz*l3*lc1*sin(r_q4) + g*lc1*m1*(-1) + g*lc1*m2*(-1) + g*l1*m4*sin(q1 - h_q2 + r_q4 + r_q5) + K_BFz*l3*lc2*sin(h_q2 - r_q4) - g*lc1*m4*sin(q1 - h_q2 + r_q4 + r_q5) + (D_BFx*l1*r_dd2*cos(2*h_q2))/2 + (D_BFx*l1*r_dd3*cos(2*h_q2))/2 - (D_BFz*l1*r_dd2*cos(2*h_q2))/2 - (D_BFz*l1*r_dd3*cos(2*h_q2))/2 - D_BFx*lc1*r_dd2*cos(2*h_q2) - D_BFx*lc1*r_dd3*cos(2*h_q2) + D_BFz*lc1*r_dd2*cos(2*h_q2) + D_BFz*lc1*r_dd3*cos(2*h_q2) + (K_BFx*l1*r_d2*cos(2*h_q2))/2 - (K_BFz*l1*r_d2*cos(2*h_q2))/2 - K_BFx*lc1*r_d2*cos(2*h_q2) + K_BFz*lc1*r_d2*cos(2*h_q2) + (K_BFx*l1*la1*sin(2*h_q2))/2 - (K_BFz*l1*la1*sin(2*h_q2))/2 + K_BFx*l1*lc1*sin(2*h_q2) - K_BFz*l1*lc1*sin(2*h_q2) - K_BFx*la1*lc1*sin(2*h_q2) + K_BFz*la1*lc1*sin(2*h_q2) + (K_BFx*l1*r_d3*sin(2*h_q2))/2 - (K_BFz*l1*r_d3*sin(2*h_q2))/2 - K_BFx*lc1*r_d3*sin(2*h_q2) + K_BFz*lc1*r_d3*sin(2*h_q2) + D_BFx*h_dq2*l1*lc1 + D_BFz*h_dq2*l1*lc1 + (D_BFx*l1*la1*r_dq4)/2 + (D_BFx*l1*la1*r_dq5)/2 + (D_BFz*l1*la1*r_dq4)/2 + (D_BFz*l1*la1*r_dq5)/2 - D_BFx*la1*lc1*r_dq4 - D_BFx*la1*lc1*r_dq5 - D_BFz*la1*lc1*r_dq4 - D_BFz*la1*lc1*r_dq5 + (K_BFx*l1*l3*sin(2*h_q2 - r_q4))/2 - (K_BFz*l1*l3*sin(2*h_q2 - r_q4))/2 - K_BFx*l3*lc1*sin(2*h_q2 - r_q4) + K_BFz*l3*lc1*sin(2*h_q2 - r_q4) + (D_BFx*l1*r_d3*r_dq4)/2 + (D_BFx*l1*r_d3*r_dq5)/2 + (D_BFz*l1*r_d3*r_dq4)/2 + (D_BFz*l1*r_d3*r_dq5)/2 - D_BFx*lc1*r_d3*r_dq4 - D_BFx*lc1*r_d3*r_dq5 - D_BFz*lc1*r_d3*r_dq4 - D_BFz*lc1*r_d3*r_dq5 - D_BFx*h_dq2*l1*lc1*cos(2*h_q2) + D_BFz*h_dq2*l1*lc1*cos(2*h_q2) + l1*l3*m4*r_ddq4*cos(h_q2 - r_q5) - l3*lc1*m4*r_ddq4*cos(h_q2 - r_q5) - (D_BFx*l1*la1*r_dq4*cos(2*h_q2))/2 - (D_BFx*l1*la1*r_dq5*cos(2*h_q2))/2 + (D_BFz*l1*la1*r_dq4*cos(2*h_q2))/2 + (D_BFz*l1*la1*r_dq5*cos(2*h_q2))/2 + D_BFx*la1*lc1*r_dq4*cos(2*h_q2) + D_BFx*la1*lc1*r_dq5*cos(2*h_q2) - D_BFz*la1*lc1*r_dq4*cos(2*h_q2) - D_BFz*la1*lc1*r_dq5*cos(2*h_q2) - (D_BFx*l1*r_d3*r_dq4*cos(2*h_q2))/2 - (D_BFx*l1*r_d3*r_dq5*cos(2*h_q2))/2 + (D_BFz*l1*r_d3*r_dq4*cos(2*h_q2))/2 + (D_BFz*l1*r_d3*r_dq5*cos(2*h_q2))/2 + D_BFx*lc1*r_d3*r_dq4*cos(2*h_q2) + D_BFx*lc1*r_d3*r_dq5*cos(2*h_q2) - D_BFz*lc1*r_d3*r_dq4*cos(2*h_q2) - D_BFz*lc1*r_d3*r_dq5*cos(2*h_q2) + (D_BFx*l1*r_d2*r_dq4*sin(2*h_q2))/2 + (D_BFx*l1*r_d2*r_dq5*sin(2*h_q2))/2 - (D_BFz*l1*r_d2*r_dq4*sin(2*h_q2))/2 - (D_BFz*l1*r_d2*r_dq5*sin(2*h_q2))/2 - D_BFx*lc1*r_d2*r_dq4*sin(2*h_q2) - D_BFx*lc1*r_d2*r_dq5*sin(2*h_q2) + D_BFz*lc1*r_d2*r_dq4*sin(2*h_q2) + D_BFz*lc1*r_d2*r_dq5*sin(2*h_q2) + (ddq1*l1*lc1*m2*cos(2*h_q2))/2 - (h_ddq2*l1*lc1*m2*cos(2*h_q2))/2 + dq1^2*l1*lc2*m2*sin(h_q2) - dq1^2*lc1*lc2*m2*sin(h_q2) - h_dq2^2*lc1*lc2*m2*sin(h_q2) - (D_BFx*l1*l3*r_dq5*cos(2*h_q2 - r_q4))/2 + (D_BFz*l1*l3*r_dq5*cos(2*h_q2 - r_q4))/2 + D_BFx*l3*lc1*r_dq5*cos(2*h_q2 - r_q4) - D_BFz*l3*lc1*r_dq5*cos(2*h_q2 - r_q4) - l1*la4*m4*r_dq5^2*sin(h_q2) + la4*lc1*m4*r_dq5^2*sin(h_q2) - l1*l3*m4*r_dq4^2*sin(h_q2 - r_q5) + l3*lc1*m4*r_dq4^2*sin(h_q2 - r_q5) + ddq1*l1*la1*m4*sin(r_q4 - h_q2 + r_q5) - ddq1*la1*lc1*m4*sin(r_q4 - h_q2 + r_q5) + D_BFz*h_dq2*l1*lc2*cos(h_q2) - D_BFz*la1*lc2*r_dq4*cos(h_q2) - D_BFz*la1*lc2*r_dq5*cos(h_q2) - D_BFz*lc2*r_d3*r_dq4*cos(h_q2) - D_BFz*lc2*r_d3*r_dq5*cos(h_q2) + (D_BFx*l1*l3*r_dq5*cos(r_q4))/2 + (D_BFz*l1*l3*r_dq5*cos(r_q4))/2 - D_BFx*l3*lc1*r_dq5*cos(r_q4) - D_BFz*l3*lc1*r_dq5*cos(r_q4) + D_BFz*lc2*r_d2*r_dq4*sin(h_q2) + D_BFz*lc2*r_d2*r_dq5*sin(h_q2) + ddq1*l1*lc2*m2*cos(h_q2) + ddq1*lc1*lc2*m2*cos(h_q2) + dq1^2*l1*la1*m4*cos(r_q4 - h_q2 + r_q5) - dq1^2*la1*lc1*m4*cos(r_q4 - h_q2 + r_q5) + h_ddq2*lc1*lc2*m2*cos(h_q2) + l1*la4*m4*r_ddq5*cos(h_q2) - la4*lc1*m4*r_ddq5*cos(h_q2) - D_BFz*l3*lc2*r_dq5*cos(h_q2 - r_q4) - 2*dq1*h_dq2*lc1*lc2*m2*sin(h_q2)
% 
%         tau2 == D_BMy*h_dq2 + I_G2z*h_ddq2 + K_BMy*h_q2 - D_BMy*r_dq4 - D_BMy*r_dq5 - K_BMy*r_q4 - K_BMy*r_q5 + h_ddq2*lc2^2*m2 - K_BFz*la4*lc2*sin(r_q4 - h_q2 + r_q5) + g*lc2*m2*sin(h_q2 + pi/2) + K_BFz*lc2*r_d3*cos(h_q2) + D_BFz*lc2*r_dd2*sin(h_q2) + D_BFz*lc2*r_dd3*sin(h_q2) - K_BFz*l1*lc2*sin(h_q2) + K_BFz*la1*lc2*sin(h_q2) + K_BFz*lc2*r_d2*sin(h_q2) + K_BFz*l3*lc2*sin(h_q2 - r_q4) + D_BFz*h_dq2*l1*lc2*cos(h_q2) - D_BFz*la1*lc2*r_dq4*cos(h_q2) - D_BFz*la1*lc2*r_dq5*cos(h_q2) - D_BFz*lc2*r_d2*r_dq4*cos(h_q2) - D_BFz*lc2*r_d2*r_dq5*cos(h_q2) + D_BFz*lc2*r_d3*r_dq4*sin(h_q2) + D_BFz*lc2*r_d3*r_dq5*sin(h_q2) - D_BFz*l3*lc2*r_dq5*cos(h_q2 - r_q4)
% 
%         tau3 == - D_BMy*h_dq2 - K_BMy*h_q2 + D_AMy*r_dq4 + D_BMy*r_dq4 + D_BMy*r_dq5 + 2*I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_AMy*r_q4 + K_BMy*r_q4 + K_BMy*r_q5 + l3^2*m4*r_ddq4 + lc3^2*m3*r_ddq4 - D_AFz*lc3*r_dd2 - K_AFz*lc3*r_d2 + la4*lc4*m4*r_ddq5 - K_BFz*lc4*r_d3*cos(h_q2) - D_BFz*lc4*r_dd2*sin(h_q2) - D_BFz*lc4*r_dd3*sin(h_q2) + K_BFz*l1*lc4*sin(h_q2) - K_BFz*la1*lc4*sin(h_q2) - K_BFz*lc4*r_d2*sin(h_q2) - K_BFx*l3*la2*sin(r_q5) - l3*m4*r_ddd3*cos(r_q4) - lc3*m3*r_ddd3*cos(r_q4) - l3*m4*r_ddd2*sin(r_q4) - lc3*m3*r_ddd2*sin(r_q4) + K_BFx*l3^2*sin(h_q2)*sin(r_q4)*sin(r_q5) - K_BFz*l3*r_d3*cos(h_q2)*cos(r_q5) + D_BFx*l3*r_dd2*cos(h_q2)*sin(r_q5) + D_BFx*l3*r_dd3*cos(h_q2)*sin(r_q5) - D_BFz*l3*r_dd2*cos(r_q5)*sin(h_q2) - D_BFz*l3*r_dd3*cos(r_q5)*sin(h_q2) - K_BFx*l1*l3*cos(h_q2)*sin(r_q5) + K_BFz*l1*l3*cos(r_q5)*sin(h_q2) + K_BFx*l3*la1*cos(h_q2)*sin(r_q5) - K_BFx*l3*la4*cos(h_q2)*sin(r_q4) + K_BFx*l3*la4*cos(r_q4)*sin(h_q2) - K_BFz*l3*la1*cos(r_q5)*sin(h_q2) + K_BFz*l3*lc4*cos(h_q2)*sin(r_q4) - K_BFz*l3*lc4*cos(r_q4)*sin(h_q2) + K_BFx*l3*r_d2*cos(h_q2)*sin(r_q5) - K_BFz*l3*r_d2*cos(r_q5)*sin(h_q2) - K_BFx*l3*r_d3*sin(h_q2)*sin(r_q5) - lc4*m4*r_ddd3*cos(r_q4)*cos(r_q5) + g*l3*m4*cos(r_q4) + g*lc3*m3*cos(r_q4) - l3*la4*m4*r_dq5^2*sin(r_q5) + l3*lc4*m4*r_dq4^2*sin(r_q5) - lc4*m4*r_ddd2*cos(r_q4)*sin(r_q5) - lc4*m4*r_ddd2*cos(r_q5)*sin(r_q4) + lc4*m4*r_ddd3*sin(r_q4)*sin(r_q5) - D_BFz*h_dq2*l1*lc4*cos(h_q2) + D_BFz*la1*lc4*r_dq4*cos(h_q2) + D_BFz*la1*lc4*r_dq5*cos(h_q2) + D_BFz*lc4*r_d2*r_dq4*cos(h_q2) + D_BFz*lc4*r_d2*r_dq5*cos(h_q2) - D_BFz*lc4*r_d3*r_dq4*sin(h_q2) - D_BFz*lc4*r_d3*r_dq5*sin(h_q2) + l3*la4*m4*r_ddq5*cos(r_q5) + l3*lc4*m4*r_ddq4*cos(r_q5) + K_BFx*l3^2*cos(h_q2)*cos(r_q4)*sin(r_q5) + K_BFz*l3^2*cos(h_q2)*cos(r_q5)*sin(r_q4) - K_BFz*l3^2*cos(r_q4)*cos(r_q5)*sin(h_q2) + K_BFz*la4*lc4*cos(h_q2)*cos(r_q4)*sin(r_q5) + K_BFz*la4*lc4*cos(h_q2)*cos(r_q5)*sin(r_q4) - K_BFz*la4*lc4*cos(r_q4)*cos(r_q5)*sin(h_q2) + K_BFz*la4*lc4*sin(h_q2)*sin(r_q4)*sin(r_q5) + g*lc4*m4*cos(r_q4)*cos(r_q5) - g*lc4*m4*sin(r_q4)*sin(r_q5) + D_BFz*l3^2*r_dq5*cos(h_q2)*cos(r_q4)*cos(r_q5) + K_BFx*l3*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFx*l3*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) + K_BFz*l3*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFz*l3*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) - D_BFx*l3^2*r_dq5*cos(h_q2)*sin(r_q4)*sin(r_q5) + D_BFx*l3^2*r_dq5*cos(r_q4)*sin(h_q2)*sin(r_q5) + D_BFz*l3^2*r_dq5*cos(r_q5)*sin(h_q2)*sin(r_q4) - D_BFz*h_dq2*l1*l3*cos(h_q2)*cos(r_q5) + D_BFz*l3*la1*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*l3*la1*r_dq5*cos(h_q2)*cos(r_q5) + D_BFz*l3*lc4*r_dq5*cos(h_q2)*cos(r_q4) + D_BFz*l3*r_d2*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*l3*r_d2*r_dq5*cos(h_q2)*cos(r_q5) - D_BFx*h_dq2*l1*l3*sin(h_q2)*sin(r_q5) + D_BFx*l3*r_d3*r_dq4*cos(h_q2)*sin(r_q5) + D_BFx*l3*r_d3*r_dq5*cos(h_q2)*sin(r_q5) - D_BFz*l3*r_d3*r_dq4*cos(r_q5)*sin(h_q2) - D_BFz*l3*r_d3*r_dq5*cos(r_q5)*sin(h_q2) + D_BFx*l3*la1*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*l3*la1*r_dq5*sin(h_q2)*sin(r_q5) + D_BFz*l3*lc4*r_dq5*sin(h_q2)*sin(r_q4) + D_BFx*l3*r_d2*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*l3*r_d2*r_dq5*sin(h_q2)*sin(r_q5) + K_BFx*l3*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFz*l3*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFx*l3*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5) + K_BFz*l3*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5)
% 
%         tau4 == - D_BMy*h_dq2 - K_BMy*h_q2 + D_BMy*r_dq4 + D_BMy*r_dq5 + I_G4z*r_ddq4 + I_G4z*r_ddq5 + K_BMy*r_q4 + K_BMy*r_q5 + K_BFz*la4*lc4*sin(r_q4 - h_q2 + r_q5) + la4*lc4*m4*r_ddq5 - lc4*m4*r_ddd3*cos(r_q4 + r_q5) - K_BFz*lc4*r_d3*cos(h_q2) - D_BFz*lc4*r_dd2*sin(h_q2) - D_BFz*lc4*r_dd3*sin(h_q2) + K_BFz*l1*lc4*sin(h_q2) - K_BFz*la1*lc4*sin(h_q2) - lc4*m4*r_ddd2*sin(r_q4 + r_q5) - K_BFz*lc4*r_d2*sin(h_q2) - K_BFz*l3*lc4*sin(h_q2 - r_q4) + g*lc4*m4*sin(pi/2 + r_q4 + r_q5) + l3*lc4*m4*r_dq4^2*sin(r_q5) - D_BFz*h_dq2*l1*lc4*cos(h_q2) + D_BFz*la1*lc4*r_dq4*cos(h_q2) + D_BFz*la1*lc4*r_dq5*cos(h_q2) + D_BFz*lc4*r_d2*r_dq4*cos(h_q2) + D_BFz*lc4*r_d2*r_dq5*cos(h_q2) - D_BFz*lc4*r_d3*r_dq4*sin(h_q2) - D_BFz*lc4*r_d3*r_dq5*sin(h_q2) + D_BFz*l3*lc4*r_dq5*cos(h_q2 - r_q4) + l3*lc4*m4*r_ddq4*cos(r_q5)
% 
%         0 == D_AFx*r_dd3 + K_AFx*r_d3 - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 + m3*r_ddd2*cos(r_q4) + m4*r_ddd2*cos(r_q4) - m3*r_ddd3*sin(r_q4) - m4*r_ddd3*sin(r_q4) - K_BFx*la2*cos(r_q5) - la4*m4*r_ddq5*sin(r_q5) + D_BFx*r_dd2*cos(h_q2)*cos(r_q5) + D_BFx*r_dd3*cos(h_q2)*cos(r_q5) - K_BFx*l1*cos(h_q2)*cos(r_q5) + K_BFx*la1*cos(h_q2)*cos(r_q5) - K_BFz*la4*cos(h_q2)*cos(r_q4) + K_BFx*r_d2*cos(h_q2)*cos(r_q5) - K_BFx*r_d3*cos(r_q5)*sin(h_q2) + K_BFz*r_d3*cos(h_q2)*sin(r_q5) + D_BFz*r_dd2*sin(h_q2)*sin(r_q5) + D_BFz*r_dd3*sin(h_q2)*sin(r_q5) - K_BFz*l1*sin(h_q2)*sin(r_q5) + K_BFz*la1*sin(h_q2)*sin(r_q5) - K_BFz*la4*sin(h_q2)*sin(r_q4) + K_BFz*r_d2*sin(h_q2)*sin(r_q5) - la4*m4*r_dq5^2*cos(r_q5) + g*m3*sin(r_q4) + g*m4*sin(r_q4) + K_BFx*la4*cos(r_q5)^2*sin(h_q2)*sin(r_q4) + K_BFz*la4*cos(r_q5)^2*sin(h_q2)*sin(r_q4) - D_BFx*h_dq2*l1*cos(r_q5)*sin(h_q2) + D_BFz*h_dq2*l1*cos(h_q2)*sin(r_q5) + D_BFx*r_d3*r_dq4*cos(h_q2)*cos(r_q5) + D_BFx*r_d3*r_dq5*cos(h_q2)*cos(r_q5) + D_BFx*la1*r_dq4*cos(r_q5)*sin(h_q2) + D_BFx*la1*r_dq5*cos(r_q5)*sin(h_q2) - D_BFz*la1*r_dq4*cos(h_q2)*sin(r_q5) - D_BFz*la1*r_dq5*cos(h_q2)*sin(r_q5) + D_BFx*r_d2*r_dq4*cos(r_q5)*sin(h_q2) + D_BFx*r_d2*r_dq5*cos(r_q5)*sin(h_q2) - D_BFz*r_d2*r_dq4*cos(h_q2)*sin(r_q5) - D_BFz*r_d2*r_dq5*cos(h_q2)*sin(r_q5) + D_BFz*r_d3*r_dq4*sin(h_q2)*sin(r_q5) + D_BFz*r_d3*r_dq5*sin(h_q2)*sin(r_q5) + K_BFx*l3*cos(h_q2)*cos(r_q4)*cos(r_q5) + K_BFx*l3*cos(r_q5)*sin(h_q2)*sin(r_q4) - K_BFz*l3*cos(h_q2)*sin(r_q4)*sin(r_q5) + K_BFz*l3*cos(r_q4)*sin(h_q2)*sin(r_q5) + K_BFx*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)^2 + K_BFz*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)^2 - D_BFx*l3*r_dq5*cos(h_q2)*cos(r_q5)*sin(r_q4) + D_BFx*l3*r_dq5*cos(r_q4)*cos(r_q5)*sin(h_q2) - D_BFz*l3*r_dq5*cos(h_q2)*cos(r_q4)*sin(r_q5) - D_BFz*l3*r_dq5*sin(h_q2)*sin(r_q4)*sin(r_q5) - K_BFx*la4*cos(h_q2)*cos(r_q5)*sin(r_q4)*sin(r_q5) + K_BFx*la4*cos(r_q4)*cos(r_q5)*sin(h_q2)*sin(r_q5) - K_BFz*la4*cos(h_q2)*cos(r_q5)*sin(r_q4)*sin(r_q5) + K_BFz*la4*cos(r_q4)*cos(r_q5)*sin(h_q2)*sin(r_q5)
%         0 == l3*m4*r_ddq4 - K_AFz*r_d2 - m3*r_ddd3*cos(r_q4) - m4*r_ddd3*cos(r_q4) - m3*r_ddd2*sin(r_q4) - m4*r_ddd2*sin(r_q4) - D_AFz*r_dd2 + lc3*m3*r_ddq4 - K_BFx*la2*sin(r_q5) + la4*m4*r_ddq5*cos(r_q5) - K_BFz*r_d3*cos(h_q2)*cos(r_q5) + D_BFx*r_dd2*cos(h_q2)*sin(r_q5) + D_BFx*r_dd3*cos(h_q2)*sin(r_q5) - D_BFz*r_dd2*cos(r_q5)*sin(h_q2) - D_BFz*r_dd3*cos(r_q5)*sin(h_q2) - K_BFx*l1*cos(h_q2)*sin(r_q5) + K_BFz*l1*cos(r_q5)*sin(h_q2) + K_BFx*la1*cos(h_q2)*sin(r_q5) - K_BFx*la4*cos(h_q2)*sin(r_q4) + K_BFx*la4*cos(r_q4)*sin(h_q2) - K_BFz*la1*cos(r_q5)*sin(h_q2) + K_BFx*r_d2*cos(h_q2)*sin(r_q5) - K_BFz*r_d2*cos(r_q5)*sin(h_q2) - K_BFx*r_d3*sin(h_q2)*sin(r_q5) + g*m3*cos(r_q4) + g*m4*cos(r_q4) - la4*m4*r_dq5^2*sin(r_q5) - D_BFz*h_dq2*l1*cos(h_q2)*cos(r_q5) + D_BFz*la1*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*la1*r_dq5*cos(h_q2)*cos(r_q5) + D_BFz*r_d2*r_dq4*cos(h_q2)*cos(r_q5) + D_BFz*r_d2*r_dq5*cos(h_q2)*cos(r_q5) - D_BFx*h_dq2*l1*sin(h_q2)*sin(r_q5) + D_BFx*r_d3*r_dq4*cos(h_q2)*sin(r_q5) + D_BFx*r_d3*r_dq5*cos(h_q2)*sin(r_q5) - D_BFz*r_d3*r_dq4*cos(r_q5)*sin(h_q2) - D_BFz*r_d3*r_dq5*cos(r_q5)*sin(h_q2) + D_BFx*la1*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*la1*r_dq5*sin(h_q2)*sin(r_q5) + D_BFx*r_d2*r_dq4*sin(h_q2)*sin(r_q5) + D_BFx*r_d2*r_dq5*sin(h_q2)*sin(r_q5) + K_BFx*l3*cos(h_q2)*cos(r_q4)*sin(r_q5) + K_BFz*l3*cos(h_q2)*cos(r_q5)*sin(r_q4) - K_BFz*l3*cos(r_q4)*cos(r_q5)*sin(h_q2) + K_BFx*l3*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFx*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFx*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) + K_BFz*la4*cos(h_q2)*cos(r_q5)^2*sin(r_q4) - K_BFz*la4*cos(r_q4)*cos(r_q5)^2*sin(h_q2) - D_BFx*l3*r_dq5*cos(h_q2)*sin(r_q4)*sin(r_q5) + D_BFx*l3*r_dq5*cos(r_q4)*sin(h_q2)*sin(r_q5) + D_BFz*l3*r_dq5*cos(r_q5)*sin(h_q2)*sin(r_q4) + K_BFx*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5) + K_BFz*la4*cos(h_q2)*cos(r_q4)*cos(r_q5)*sin(r_q5) + K_BFx*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + K_BFz*la4*cos(r_q5)*sin(h_q2)*sin(r_q4)*sin(r_q5) + D_BFz*l3*r_dq5*cos(h_q2)*cos(r_q4)*cos(r_q5)
% 
%         ];

% Solve in steps
% eqns = [%tau_intAz - tau_intBz + I_G1z*ddq1 - lc2*(f_intBj + m2gj - m2*(l1*sin(h_q2)*dq1^2 + lc2*(ddq1 + h_ddq2) + ddq1*l1*cos(h_q2))) - (l1 - lc1)*(sin(h_q2)*(m4gi - f_intBi + m4*(la4*r_dq5^2 + cos(r_q5)*(l3*r_dq4^2 + sin(r_q4)*(r_ddd3 + ddq1*la1) - cos(r_q4)*(r_ddd2 - dq1^2*la1)) + sin(r_q5)*(cos(r_q4)*(r_ddd3 + ddq1*la1) - l3*r_ddq4 + sin(r_q4)*(r_ddd2 - dq1^2*la1)))) - cos(h_q2)*(f_intBj - m4gj + m4*(sin(r_q5)*(l3*r_dq4^2 + sin(r_q4)*(r_ddd3 + ddq1*la1) - cos(r_q4)*(r_ddd2 - dq1^2*la1)) + la4*r_ddq5 - cos(r_q5)*(cos(r_q4)*(r_ddd3 + ddq1*la1) - l3*r_ddq4 + sin(r_q4)*(- la1*dq1^2 + r_ddd2))))) + I_G2z*(ddq1 + h_ddq2) - lc1*(f_intAj + m1gj + sin(h_q2)*(f_intBi + m2gi + m2*(lc2*(dq1 + h_dq2)^2 - h_ddq2*l1*sin(h_q2) + dq1^2*l1*cos(h_q2))) + cos(h_q2)*(f_intBj + m2gj - m2*(l1*sin(h_q2)*dq1^2 + lc2*(ddq1 + h_ddq2) + ddq1*l1*cos(h_q2))) - ddq1*lc1*m1)
% 
%         tau2 == I_G2z*(ddq1 + h_ddq2) - lc2*(f_intBj + m2gj - m2*(l1*sin(h_q2)*dq1^2 + lc2*(ddq1 + h_ddq2) + ddq1*l1*cos(h_q2))) - tau_intBz
% 
%         tau3 == tau_intAz + tau_intBz + 2*I_G4z*ddq1 + 2*I_G4z*r_ddq4 + I_G4z*r_ddq5 + f_intAj*lc3 + f_intBj*lc4 - lc3*m3gj - lc4*m4gj + f_intBj*l3*cos(r_q5) + l3^2*m4*r_ddq4 + lc3^2*m3*r_ddq4 - l3*m4gj*cos(r_q5) + f_intBi*l3*sin(r_q5) - l3*m4gi*sin(r_q5) + la4*lc4*m4*r_ddq5 - lc4*m4*r_ddd3*cos(r_q4 + r_q5) - lc4*m4*r_ddd2*sin(r_q4 + r_q5) - l3*m4*r_ddd3*cos(r_q4) - lc3*m3*r_ddd3*cos(r_q4) - l3*m4*r_ddd2*sin(r_q4) - lc3*m3*r_ddd2*sin(r_q4) + dq1^2*la1*lc4*m4*sin(r_q4 + r_q5) + dq1^2*l3*la1*m4*sin(r_q4) + dq1^2*la1*lc3*m3*sin(r_q4) - l3*la4*m4*r_dq5^2*sin(r_q5) + l3*lc4*m4*r_dq4^2*sin(r_q5) - ddq1*la1*lc4*m4*cos(r_q4 + r_q5) - ddq1*l3*la1*m4*cos(r_q4) - ddq1*la1*lc3*m3*cos(r_q4) + l3*la4*m4*r_ddq5*cos(r_q5) + l3*lc4*m4*r_ddq4*cos(r_q5)
% 
%         tau4 == la1*lc4*m4*sin(r_q4 + r_q5)*dq1^2 + l3*lc4*m4*sin(r_q5)*r_dq4^2 + tau_intBz + I_G4z*ddq1 + I_G4z*r_ddq4 + I_G4z*r_ddq5 + f_intBj*lc4 - lc4*m4gj + la4*lc4*m4*r_ddq5 - lc4*m4*r_ddd3*cos(r_q4 + r_q5) - lc4*m4*r_ddd2*sin(r_q4 + r_q5) - ddq1*la1*lc4*m4*cos(r_q4 + r_q5) + l3*lc4*m4*r_ddq4*cos(r_q5)
% 
%         0 == f_intAi - m3gi + f_intBi*cos(r_q5) - m4gi*cos(r_q5) - f_intBj*sin(r_q5) + m4gj*sin(r_q5) - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 + m3*r_ddd2*cos(r_q4) + m4*r_ddd2*cos(r_q4) - m3*r_ddd3*sin(r_q4) - m4*r_ddd3*sin(r_q4) - ddq1*la1*m3*sin(r_q4) - ddq1*la1*m4*sin(r_q4) - la4*m4*r_ddq5*sin(r_q5) - dq1^2*la1*m3*cos(r_q4) - dq1^2*la1*m4*cos(r_q4) - la4*m4*r_dq5^2*cos(r_q5)
%         0 == f_intAj - m3gj + f_intBj*cos(r_q5) - m4gj*cos(r_q5) + f_intBi*sin(r_q5) - m4gi*sin(r_q5) - m3*r_ddd3*cos(r_q4) - m4*r_ddd3*cos(r_q4) - m3*r_ddd2*sin(r_q4) - m4*r_ddd2*sin(r_q4) + l3*m4*r_ddq4 + lc3*m3*r_ddq4 - ddq1*la1*m3*cos(r_q4) - ddq1*la1*m4*cos(r_q4) + la4*m4*r_ddq5*cos(r_q5) + dq1^2*la1*m3*sin(r_q4) + dq1^2*la1*m4*sin(r_q4) - la4*m4*r_dq5^2*sin(r_q5)
% 
%         ];

% Solve in steps with items collected
% eqns = [
%         tau2 == (m2*lc2^2 + I_G2z)*h_ddq2 - tau_intBz - lc2*(f_intBj + m2gj)
% 
%         tau3 == l3*lc4*m4*sin(r_q5)*r_dq4^2 - l3*la4*m4*sin(r_q5)*r_dq5^2 + tau_intAz + tau_intBz - r_ddd3*(l3*m4*cos(r_q4) + lc3*m3*cos(r_q4) + lc4*m4*cos(r_q4 + r_q5)) + f_intAj*lc3 + f_intBj*lc4 - lc3*m3gj - lc4*m4gj - 
%                   r_ddd2*(l3*m4*sin(r_q4) + lc3*m3*sin(r_q4) + lc4*m4*sin(r_q4 + r_q5)) + r_ddq4*(m4*l3^2 + lc4*m4*cos(r_q5)*l3 + m3*lc3^2 + 2*I_G4z) + r_ddq5*(I_G4z + la4*lc4*m4 + l3*la4*m4*cos(r_q5)) + 
%                   f_intBj*l3*cos(r_q5) - l3*m4gj*cos(r_q5) + f_intBi*l3*sin(r_q5) - l3*m4gi*sin(r_q5)
% 
%         tau4 == l3*lc4*m4*sin(r_q5)*r_dq4^2 + tau_intBz + r_ddq5*(I_G4z + la4*lc4*m4) + f_intBj*lc4 - lc4*m4gj + r_ddq4*(I_G4z + l3*lc4*m4*cos(r_q5)) - lc4*m4*r_ddd3*cos(r_q4 + r_q5) - lc4*m4*r_ddd2*sin(r_q4 + r_q5)
% 
%         0 == (m3*cos(r_q4) + m4*cos(r_q4))*r_ddd2 + (- m3*sin(r_q4) - m4*sin(r_q4))*r_ddd3 + (-la4*m4*sin(r_q5))*r_ddq5 + f_intAi - m3gi + f_intBi*cos(r_q5) - m4gi*cos(r_q5) - f_intBj*sin(r_q5) + m4gj*sin(r_q5) - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 - la4*m4*r_dq5^2*cos(r_q5)
%         0 == - la4*m4*sin(r_q5)*r_dq5^2 + f_intAj - m3gj - r_ddd3*(m3*cos(r_q4) + m4*cos(r_q4)) - r_ddd2*(m3*sin(r_q4) + m4*sin(r_q4)) + f_intBj*cos(r_q5) - m4gj*cos(r_q5) + r_ddq4*(l3*m4 + lc3*m3) + f_intBi*sin(r_q5) - m4gi*sin(r_q5) + la4*m4*r_ddq5*cos(r_q5)
%         ];


% test
syms a1 a2
syms b1 b2 b3 b4 b5
syms c1 c2 c3 c4 c5
syms d1 d2 d3 d5
syms e1 e2 e3 e4 e5
    
eqns = [
    tau2 == a1 + a2*h_ddq2

    0 == b1 + (b2)*r_ddd2 + (b3)*r_ddd3 + (b4)*r_ddq4 + (b5)*r_ddq5  % no active input

    tau4 == c1 + (c2)*r_ddd2 + (c3)*r_ddd3 + (c4)*r_ddq4 + (c5)*r_ddq5

    0 == d1 + (d2)*r_ddd2 + (d3)*r_ddd3 + (d5)*r_ddq5
    0 == e1 + (e2)*r_ddd2 + (e3)*r_ddd3 + (e4)*r_ddq4 + (e5)*r_ddq5
    ];


%tobeElimitedVars = [];
%expr = eliminate(eqns, tobeElimitedVars);
ddthetaSolved = solve(eqns, [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5]);
%ddthetaSolved = solve(expr, [ddq_2, ddq_3, ddq_4]);


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
vSym.l1 = 0.4; vSym.l2 = 0.4; vSym.l3 = 0.2; vSym.l4 = 0.2;
vSym.lc1 = 0.5; vSym.la1 = 0.2; vSym.lb1 = 0.2; vSym.lc2 = 0.5; 
vSym.lc3 = 0.173; vSym.lc4 = 0.1; vSym.la4 = 0.2; vSym.la2 = 0.2; 
vSym.K_AFz = 3000; vSym.K_AFx = -6000; vSym.K_AMy = 50;
vSym.K_BFz = 3000; vSym.K_BFx = 6000; vSym.K_BMy = 50;
vSym.D_AFz = 150; vSym.D_AFx = -150; vSym.D_AMy = 15;
vSym.D_BFz = 150; vSym.D_BFx = 150; vSym.D_BMy = 15;
% vSym.K_AFz = 0; vSym.K_AFx = 0; vSym.K_AMy = 0;
% vSym.K_BFz = 0; vSym.K_BFx = 0; vSym.K_BMy = 0;
% vSym.D_AFz = 1000; vSym.D_AFx = -1000; vSym.D_AMy = 10;
% vSym.D_BFz = 0; vSym.D_BFx = 0; vSym.D_BMy = 0;
% vSym.K_AMy = 0; vSym.K_BMy = 0;
% vSym.D_AFz = 300; vSym.D_AFx = 300; vSym.D_AMy = 20;
% vSym.D_BFz = 300; vSym.D_BFx = 300; vSym.D_BMy = 20;

%% solve the dynamic equations
tspan = [0 4];
%tspan = 0:0.05:1;
y0 = [q1, dq1, -pi/2, 0, 0, 0, 0, 0, 0, 0, -pi/2, 0];  %initial condition
tau = [0, 0, 0, 0];
[t,y] = ode23(@(t,y) runrobot(t,y, tau, vSym), tspan, y0);

%% Plot trajectory
% figure('Renderer', 'painters', 'Position', [300 300 800 800])
% plot(t, y(:, [3, 11]))
% legend(["h th", "r th"])

%% Visualise
dt = 0.02;
q1 = y(:, 1); h_q2 = y(:, 3); r_d2 = y(:, 5); r_d3 = y(:, 7); r_q4 = y(:, 9); r_q5 = y(:, 11);
dq1 = y(:, 2); h_dq2 = y(:, 4); r_dd2 = y(:, 6); r_dd3 = y(:, 8); r_dq4 = y(:, 10); r_dq5 = y(:, 12);

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
    

figure('Renderer', 'painters', 'Position', [300 300 800 700])
for i = 1:length(y)
    [j1, h_j2, h_j3, r_j2, r_j3, r_j4, r_j5, r_j6, T] = calc_joint_position(vSym, q1(i), h_q2(i), r_d2(i), r_d3(i), r_q4(i), r_q5(i));

    clf
    plot_3d_points({j1, h_j2, h_j3});
    plot_3d_points({j1, r_j2, r_j3, r_j4, r_j5, r_j6});
    view(0, 90)

    axis([-0.4,1,-0.8,0.6, -2,4])
    pause(dt)
end

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
        vVar.tau2 = 0;
        vVar.tau4 = 5*sin(2*pi*t);
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

function [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5] = calc_equation_of_motion(vSym, vVar)
    % syms l1 l2 l3 l4
    % syms lc1 la1 lb1 lc2  lc3 lc4 la4 la2
    % syms q1 h_q2 dq1 h_dq2   
    % syms r_d2 r_dd2 
    % syms r_d3 r_dd3 
    % syms r_q4 r_dq4 
    % syms r_q5 r_dq5 
    % syms K_AFz K_AFx K_AMy
    % syms K_BFz K_BFx K_BMy
    % syms D_AFz D_AFx D_AMy
    % syms D_BFz D_BFx D_BMy
    % syms m1 m2 m3 m4 g
    % syms I_G1z I_G2z I_G3z I_G4z
    % syms tau1 tau2 tau3 tau4
    % symbols = {m1, m2, m3, m4, g, I_G1z, I_G2z, I_G3z, I_G4z, ...
    %        l1, l2, l3, l4,...
    %        lc1, la1, lb1, lc2, lc3, lc4, la4, la2,...
    %        K_AFz, K_AFx, K_AMy,...
    %        K_BFz, K_BFx, K_BMy, ...
    %        D_AFz, D_AFx, D_AMy,...
    %        D_BFz, D_BFx, D_BMy, ...
    %        q1, dq1, ...
    %        h_q2, h_dq2,...
    %        r_d2, r_dd2,...
    %        r_d3, r_dd3,...
    %        r_q4, r_dq4,...
    %        r_q5, r_dq5,...
    %        tau1, tau2, tau3, tau4};
    % symbol_vals = {vSym.m1, vSym.m2, vSym.m3, vSym.m4, vSym.g,...
    %                vSym.I_G1z, vSym.I_G2z, vSym.I_G3z, vSym.I_G4z, ...
    %                vSym.l1, vSym.l2, vSym.l3, vSym.l4,...
    %                vSym.lc1, vSym.la1, vSym.lb1, vSym.lc2, ...
    %                vSym.lc3, vSym.lc4, vSym.la4, vSym.la2, ...
    %                vSym.K_AFz, vSym.K_AFx, vSym.K_AMy,...
    %                vSym.K_BFz, vSym.K_BFx, vSym.K_BMy, ...
    %                vSym.D_AFz, vSym.D_AFx, vSym.D_AMy,...
    %                vSym.D_BFz, vSym.D_BFx, vSym.D_BMy};
    % variable_vals = {vVar.q1, vVar.dq1, ...
    %                  vVar.h_q2, vVar.h_dq2,...
    %                  vVar.r_d2, vVar.r_dd2,...
    %                  vVar.r_d3, vVar.r_dd3,...
    %                  vVar.r_q4, vVar.r_dq4,...
    %                  vVar.r_q5, vVar.r_dq5,...
    %                  vVar.tau1, vVar.tau2, vVar.tau3, vVar.tau4};

    % Assign values to symbols
    l1 = vSym.l1;  l2 = vSym.l2;  l3 = vSym.l3;  l4 = vSym.l4;
    lc1 = vSym.lc1;  la1 = vSym.la1;  lb1 = vSym.lb1;  lc2 = vSym.lc2;
    lc3 = vSym.lc3;  lc4 = vSym.lc4;  la4 = vSym.la4;  la2 = vSym.la2;
    K_AFz = vSym.K_AFz;  K_AFx = vSym.K_AFx;  K_AMy = vSym.K_AMy;
    K_BFz = vSym.K_BFz;  K_BFx = vSym.K_BFx;  K_BMy = vSym.K_BMy;
    D_AFz = vSym.D_AFz;  D_AFx = vSym.D_AFx;  D_AMy = vSym.D_AMy;
    D_BFz = vSym.D_BFz;  D_BFx = vSym.D_BFx;  D_BMy = vSym.D_BMy;
    m1 = vSym.m1;  m2 = vSym.m2;  m3 = vSym.m3;  m4 = vSym.m4;  g = vSym.g;
    I_G1z = vSym.I_G1z;  I_G2z = vSym.I_G2z;  I_G3z = vSym.I_G3z;  I_G4z = vSym.I_G4z;

    q1 = vVar.q1;  h_q2 = vVar.h_q2;  dq1 = vVar.dq1;  h_dq2 = vVar.h_dq2;
    r_d2 = vVar.r_d2;  r_dd2 = vVar.r_dd2;
    r_d3 = vVar.r_d3;  r_dd3 = vVar.r_dd3;
    r_q4 = vVar.r_q4;  r_dq4 = vVar.r_dq4;
    r_q5 = vVar.r_q5;  r_dq5 = vVar.r_dq5;
    tau1 = vVar.tau1; tau2 = vVar.tau2; tau3 = vVar.tau3; tau4 = vVar.tau4; 

    f_intBi = D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))));
    f_intBj = - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))));
    m4gi = m4*g*cos(q1+r_q4+r_q5);
    m4gj = -m4*g*sin(q1+r_q4+r_q5);
    tau_intBz  = K_BMy*(r_q4 + r_q5 - h_q2) + D_BMy*(r_dq4 + r_dq5 - h_dq2);
    f_intAi  = - K_AFx * r_d2 - D_AFx * r_dd2;
    f_intAj = - K_AFz * r_d3 - D_AFz * r_dd3;
    m3gi  = m3*g*cos(q1+r_q4);
    m3gj  = -m3*g*sin(q1+r_q4);
    tau_intAz  = K_AMy*r_q4 + D_AMy*r_dq4;
    %m2gi  = subs(m2*g*cos(q1+h_q2), symbols, [symbol_vals, variable_vals]);
    m2gj  = -m2*g*sin(q1+h_q2);
    %m1gi  = subs(m1*g*cos(q1), symbols, [symbol_vals, variable_vals]);
    %m1gj   = subs(-m1*g*sin(q1), symbols, [symbol_vals, variable_vals]);

    % f_intBi = - f_intBi;
    % f_intBj = -f_intBj;
    % f_intAi = - f_intAi;
    % f_intAj = - f_intAj;
    % tau_intAz = - tau_intAz;
    % tau_intBz = - tau_intBz;

    a1 = - tau_intBz - lc2*(f_intBj + m2gj);
    a2 = m2*lc2^2 + I_G2z;
    b1 = l3*lc4*m4*sin(r_q5)*r_dq4^2 - l3*la4*m4*sin(r_q5)*r_dq5^2 + tau_intAz + tau_intBz + f_intAj*lc3 + f_intBj*lc4 - lc3*m3gj - lc4*m4gj + f_intBj*l3*cos(r_q5) - l3*m4gj*cos(r_q5) + f_intBi*l3*sin(r_q5) - l3*m4gi*sin(r_q5);
    b2 = - (l3*m4*sin(r_q4) + lc3*m3*sin(r_q4) + lc4*m4*sin(r_q4 + r_q5));
    b3 = - (l3*m4*cos(r_q4) + lc3*m3*cos(r_q4) + lc4*m4*cos(r_q4 + r_q5));
    b4 = m4*l3^2 + lc4*m4*cos(r_q5)*l3 + m3*lc3^2 + 2*I_G4z;
    b5 = I_G4z + la4*lc4*m4 + l3*la4*m4*cos(r_q5);
    c1 = l3*lc4*m4*sin(r_q5)*r_dq4^2 + tau_intBz + f_intBj*lc4 - lc4*m4gj;
    c2 = - lc4*m4*sin(r_q4 + r_q5);
    c3 = - lc4*m4*cos(r_q4 + r_q5);
    c4 = I_G4z + l3*lc4*m4*cos(r_q5);
    c5 = I_G4z + la4*lc4*m4;
    d1 =  + f_intAi - m3gi + f_intBi*cos(r_q5) - m4gi*cos(r_q5) - f_intBj*sin(r_q5) + m4gj*sin(r_q5) - l3*m4*r_dq4^2 - lc3*m3*r_dq4^2 - la4*m4*r_dq5^2*cos(r_q5);
    d2 = m3*cos(r_q4) + m4*cos(r_q4);
    d3 = - m3*sin(r_q4) - m4*sin(r_q4);
    d5 = -la4*m4*sin(r_q5);
    e1 = - la4*m4*sin(r_q5)*r_dq5^2 + f_intAj - m3gj + f_intBj*cos(r_q5) - m4gj*cos(r_q5) + f_intBi*sin(r_q5) - m4gi*sin(r_q5);
    e2 = -(m3*sin(r_q4) + m4*sin(r_q4));
    e3 = -(m3*cos(r_q4) + m4*cos(r_q4));
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


function [j1, h_j2, h_j3, r_j2, r_j3, r_j4, r_j5, r_j6, T] = calc_joint_position(vSym, q1, h_q2, r_d2, r_d3, r_q4, r_q5)
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
    scatter3(xs, ys, zs, 'o'); % 'o' is for circle marker

    % Hold on to add lines to the same plot
    hold on;

    % Connect the points with lines
    plot3(xs, ys, zs);

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
