
syms l1 l2 l3 l4
syms lc1 la1 lc2  lc3 lc4 la4 la2
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

% Matrix multiplicationd
H_T01 = calc_homotransmtx(0, 0, 0, -pi/2 + q1);
H_T12 = calc_homotransmtx(l1, 0, 0, h_q2);
H_T23a = calc_homotransmtx(la2, 0, 0, 0);
H_T02 = H_T01 * H_T12;
H_T03a = H_T02 * H_T23a;

R_T01 = calc_homotransmtx(0, 0, 0, q1);
R_T12 = calc_homotransmtx(0, pi/2, la1 + r_d2, pi/2);
R_T23 = calc_homotransmtx(0, -pi/2, r_d3, -pi/2);
R_T34 = calc_homotransmtx(0, -pi/2, 0, r_q4);
R_T45 = calc_homotransmtx(l3, 0, 0, r_q5);
R_T56a = calc_homotransmtx(la4, 0, 0, 0);
R_T02 = R_T01 * R_T12;
R_T03 = R_T02 * R_T23;
R_T04 = R_T03 * R_T34;
R_T05 = R_T04 * R_T45;
R_T06a = R_T05 * R_T56a;


ae2 = [0; la1*ddq1; -la1*dq1^2 + r_ddd2];
R23 = get_rotmax(R_T23);
R32 = R23.';
re3 = [0; 0; 0];
dw3 = [0; 0; 0];
w3 = [0; 0; 0];
ae3 = R32*ae2 + cross(dw3, re3) + cross(w3, cross(w3, re3)) + [0; 0; r_ddd3];

R34 = get_rotmax(R_T34);
R43 = R34.';
re4 = [l3; 0; 0];
rc4 = [lc3; 0; 0];
dw4 = [0; 0; r_ddq4];
w4 = [0; 0; r_dq4];
ae4 = R43*ae3 + cross(dw4, re4) + cross(w4, cross(w4, re4));
ac4 = R43*ae3 + cross(dw4, rc4) + cross(w4, cross(w4, rc4));

R45 = get_rotmax(R_T45);
R54 = R45.';
rc5 = [la4; 0; 0];
dw5 = [0; 0; r_ddq5];
w5 = [0; 0; r_dq5];
ac5 = R54*ae4 + cross(dw5, rc5) + cross(w5, cross(w5, rc5));

% interaction force
% position of attachment point B in frame {2}
H3aR6a_T = inv_homotransmtx(H_T03a) * R_T06a;

H_R01 = get_rotmax(H_T01);
H_R02 = get_rotmax(H_T02);
h_z01 = H_R01 * [0; 0; 1];
h_z02 = H_R02 * [0; 0; 1];
h_o1 = H_T01(1:3, 4);
h_o2 = H_T02(1:3, 4);
H_Jv02 = [cross(h_z01, h_o1), cross(h_z02, h_o2)];
H_Jw02 = [h_z01, h_z02];
H_J02 = [H_Jv02; H_Jw02];
vLink0_2 = H_J02 * [dq1; h_dq2];
vLink2_2_v = H_R02.' * vLink0_2(1:3, 1);
vLink2_2_w = H_R02.' * vLink0_2(4:6, 1);
vLink2_2 = simplify([vLink2_2_v; vLink2_2_w]);

R_R01 = get_rotmax(R_T01);
R_R02 = get_rotmax(R_T02);
R_R03 = get_rotmax(R_T03);
R_R04 = get_rotmax(R_T04);
R_R05 = get_rotmax(R_T05);
r_z01 = R_R01 * [0; 0; 1];
r_z02 = R_R02 * [0; 0; 1];
r_z03 = R_R03 * [0; 0; 1];
r_z04 = R_R04 * [0; 0; 1];
r_z05 = R_R05 * [0; 0; 1];
r_o1 = R_T01(1:3, 4);
r_o2 = R_T02(1:3, 4);
r_o3 = R_T03(1:3, 4);
r_o4 = R_T04(1:3, 4);
r_o5 = R_T05(1:3, 4);
R_Jv05 = [cross(r_z01, r_o1), r_z02, r_z02, cross(r_z04, r_o4), cross(r_z05, r_o5)];
R_Jw05 = [r_z01, [0;0;0], [0;0;0], r_z04, r_z05];
R_J05 = [R_Jv05; R_Jw05];
vLink0_4 = R_J05 * [dq1; r_dd2; r_dd3; r_dq4; r_dq5];
vLink2_4_v = H_R02.' * vLink0_4(1:3, 1);
vLink2_4_w = H_R02.' * vLink0_4(4:6, 1);
vLink2_4 = [vLink2_4_v; vLink2_4_w];
vLink2_link4_inFrame2 = simplify(vLink2_4 - vLink2_2);

%% 
% f_intB = H3aR6a_T(1:3, 4) .* [K_BFx; K_BFz; 0] + vLink2_link4_inFrame2(1:3, 1) .* [D_BFx; D_BFz; 0];
syms f_intBi  % D_BFx*(r_dd2*cos(h_q2) + r_dd3*cos(h_q2) - h_dq2*l1*sin(h_q2) + r_d3*r_dq4*cos(h_q2) + r_d3*r_dq5*cos(h_q2) + la1*r_dq4*sin(h_q2) + la1*r_dq5*sin(h_q2) + r_d2*r_dq4*sin(h_q2) + r_d2*r_dq5*sin(h_q2) + l3*r_dq5*sin(h_q2 - r_q4)) - K_BFx*((la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + (la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))) - (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))))
syms f_intBj  % - D_BFz*(r_dd2*sin(h_q2) + r_dd3*sin(h_q2) - la1*r_dq4*cos(h_q2) - la1*r_dq5*cos(h_q2) - r_d2*r_dq4*cos(h_q2) - r_d2*r_dq5*cos(h_q2) + r_d3*r_dq4*sin(h_q2) + r_d3*r_dq5*sin(h_q2) - l3*r_dq5*cos(h_q2 - r_q4) + h_dq2*l1*cos(h_q2)) - K_BFz*((la2*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + l1*sin(q1 - pi/2))*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) - (la2*(cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2)) + l1*cos(q1 - pi/2))*(cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2)) + (cos(h_q2)*sin(q1 - pi/2) + sin(h_q2)*cos(q1 - pi/2))*(l3*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(q1)*(la1 + r_d2) - r_d3*cos(q1) + la4*(cos(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)) + sin(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)))) + (cos(h_q2)*cos(q1 - pi/2) - sin(h_q2)*sin(q1 - pi/2))*(l3*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) + cos(q1)*(la1 + r_d2) + r_d3*sin(q1) + la4*(cos(r_q5)*(cos(q1)*cos(r_q4) - sin(q1)*sin(r_q4)) - sin(r_q5)*(cos(q1)*sin(r_q4) + cos(r_q4)*sin(q1)))))
f_intB = [f_intBi; f_intBj; 0];

%m4g = [m4*g*cos(q1+r_q4+r_q5); -m4*g*sin(q1+r_q4+r_q5); 0];
syms m4gi  % m4*g*cos(q1+r_q4+r_q5)
syms m4gj  % -m4*g*sin(q1+r_q4+r_q5)
m4g = [m4gi; m4gj; 0];

f4 = simplify(m4*ac5 - m4g + f_intB);
dhc4 = [0; 0; I_G4z*(ddq1 + r_ddq4 + r_ddq5)]; % rate of change of angular momentum of body 4

% tau_intB = [0; 0; K_BMy*(r_q4 + r_q5 - h_q2) + D_BMy*(r_dq4 + r_dq5 - h_dq2)];
syms tau_intBz  % K_BMy*(r_q4 + r_q5 - h_q2) + D_BMy*(r_dq4 + r_dq5 - h_dq2)
tau_intB = [0; 0; tau_intBz];

tau4 = simplify(dhc4 - cross(f4, [lc4;0;0]) + tau_intB);

%% 
% f_intA = [K_AFx * r_d2 + D_AFx * r_dd2; + K_AFz * r_d3 + D_AFz * r_dd3; 0];
syms f_intAi  % K_AFx * r_d3 + D_AFx * r_dd3
syms f_intAj  % + K_AFz * r_d2 + D_AFz * r_dd2
f_intA = [f_intAi; f_intAj; 0];

% m3g = [m3*g*cos(q1+r_q4); -m3*g*sin(q1+r_q4); 0];
syms m3gi  % m3*g*cos(q1+r_q4)
syms m3gj  % -m3*g*sin(q1+r_q4)
m3g = [m3gi; m3gj; 0];

f3 = simplify(m3 * ac4 + R45*f4 - m3g + f_intA);  % equals zero

% tau_intA = [0; 0; K_AMy*r_q4 + D_AMy*r_dq4];
syms tau_intAz  % K_AMy*r_q4 + D_AMy*r_dq4
tau_intA = [0; 0; tau_intAz];

dhc3 = [0; 0; I_G4z*(ddq1 + r_ddq4)]; % rate of change of angular momentum of body 3
tau3 = simplify(R45*tau4 - cross(f3, [lc3;0;0]) - cross(R45*f4, [l3-lc3;0;0]) + dhc3 + tau_intA);  % equals zero, R45 is frame 45, means rigid body 34


%% 
% Human side
h_ac2 = [-lc2*(dq1 + h_dq2)^2 - l1*cos(h_q2)*dq1^2  + l1*sin(h_q2)*h_ddq2;...
         lc2*(ddq1 + h_ddq2) + l1*sin(h_q2)*dq1^2 + l1*cos(h_q2)*ddq1;...
         0];

% m2g = [m2*g*cos(q1+h_q2); -m2*g*sin(q1+h_q2); 0];
syms m2gi  % m2*g*cos(q1+h_q2)
syms m2gj  % -m2*g*sin(q1+h_q2)
m2g = [m2gi; m2gj; 0];

f2 = simplify(m2*h_ac2 - m2g - f_intB);

dhc2 = [0; 0; I_G2z*(ddq1 + h_ddq2)]; % rate of change of angular momentum of body 2
tau2 = simplify(dhc2 - cross(f2, [lc2;0;0]) - tau_intB);


%% 
h_ac1 = [-lc1*dq1^2; lc1*ddq1; 0];

% m1g = [m1*g*cos(q1); -m1*g*sin(q1); 0];
syms m1gi  % m1*g*cos(q1)
syms m1gj  % -m1*g*sin(q1)
m1g = [m1gi; m1gj; 0];

R12 = get_rotmax(H_T12);
f1 = m1*h_ac1 + R12*f2 - m1g - f_intA;

dhc1 = [0; 0; I_G1z*ddq1]; % rate of change of angular momentum of body 1
tau1 = simplify(R12*tau2 - cross(f1, [lc1;0;0]) - cross(R12*f2, [l1-lc1;0;0]) + dhc1 - tau_intA);


%% Collect terms
% my_variables = [q1, dq1, ddq1, ...
%                 h_q2, h_dq2, h_ddq2, ...
%                 r_d2, r_dd2, r_ddd2, ...
%                 r_d3, r_dd3, r_ddd3, ...
%                 r_q4, r_dq4, r_ddq4, ...
%                 r_q5, r_dq5, r_ddq5];
my_variables = [h_ddq2, r_ddd2, r_ddd3, r_ddq4, r_ddq5];


tau1_eqn = subs(tau1, {q1, dq1, ddq1}, {pi/2, 0, 0});
tau1_eqn = collect(tau1_eqn, my_variables);

tau2_eqn = subs(tau2, {q1, dq1, ddq1}, {pi/2, 0, 0});
tau2_eqn = collect(tau2_eqn, my_variables)

tau3_eqn = subs(tau3, {q1, dq1, ddq1}, {pi/2, 0, 0});
tau3_eqn = collect(tau3_eqn, my_variables)

tau4_eqn = subs(tau4, {q1, dq1, ddq1}, {pi/2, 0, 0});
tau4_eqn = collect(tau4_eqn, my_variables)

f3x_eqn = subs(f3(1), {q1, dq1, ddq1}, {pi/2, 0, 0});
f3x_eqn = collect(f3x_eqn, my_variables)

f3y_eqn = subs(f3(2), {q1, dq1, ddq1}, {pi/2, 0, 0});
f3y_eqn = collect(f3y_eqn, my_variables)

% tau2_eqn_chr = latex(tau2_eqn(3))
% tau3_eqn_chr = latex(tau3_eqn(3))
% tau4_eqn_chr = latex(tau4_eqn(3))
% f3x_eqn_chr = latex(f3x_eqn)
% f3y_eqn_chr = latex(f3x_eqn)


%% Calculate Tl1
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
        if cos_alpha < 1e-10
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
        if cos_theta < 1e-10
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

function [inv_T] = inv_homotransmtx(T)
    R = T(1:3, 1:3);
    d = T(1:3, 4);
    inv_T = [R.', -R.'*d;
             0,0,0,1];
end

function [R] = get_rotmax(T)
    R = T(1:3, 1:3);
end
