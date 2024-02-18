
syms l1 l2 l3 l4
syms lc1 la1 lb1 lc2  lc3 lc4 la4 la2
syms q1 h_q2 dq1 h_dq2 ddq1 h_ddq2 
syms r_d2 r_dd2 r_ddd2
syms r_d3 r_dd3 r_ddd3
syms r_q4 r_dq4 r_ddq4
syms r_q5 r_dq5 r_ddq5
syms K_AFz K_AFx K_AMy
syms K_BFz K_BFx K_BMy
syms m1 m2 m3 m4 g
syms I_G1z I_G2z I_G3z I_G4z

% Matrix multiplication
t_int = ks * (h_q2 - r_q5);
h_ddq2 = g / (3*l2) * sin(h_q2) + t_int;
r_ddq5 = g / (3*l4) * sin(r_q5) - t_int;

