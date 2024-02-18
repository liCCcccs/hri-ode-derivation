
syms r_l1 r_l2 
syms r_lc1 r_lc2 
syms h_l1 h_l2
syms h_lc1 h_lc2

syms r_q1 r_dq1 r_ddq1
syms r_q2 r_dq2 r_ddq2
syms h_q1 h_dq1 h_ddq1 
syms h_q2 h_dq2 h_ddq2 

syms r_m1 r_m2 
syms h_m1 h_m2 g

syms r_I_G1z r_I_G2z
syms h_I_G1z h_I_G2z

syms Ks Ds As

% Matrix multiplication
H_T01 = calc_homotransmtx(0, 0, 0, -pi/2 + h_q1);
H_T12 = calc_homotransmtx(h_l1, 0, 0, h_q2);
H_T02 = H_T01 * H_T12;

R_T01 = calc_homotransmtx(0, 0, 0, -pi/2 + r_q1);
R_T12 = calc_homotransmtx(r_l1, 0, 0, r_q2);
R_T02 = R_T01 * R_T12;

%%
tau_int = [0; 0; Ks * (r_q2 - h_q2) + Ds * (r_dq2 - h_dq2) + As * r_q2];

%% 
% Robot side
r_h_ac2 = [-r_lc2*(r_dq1 + r_dq2)^2 - r_l1*cos(r_q2)*r_dq1^2  + r_l1*sin(r_q2)*r_ddq2;...
         r_lc2*(r_ddq1 + r_ddq2) + r_l1*sin(r_q2)*r_dq1^2 + r_l1*cos(r_q2)*r_ddq1;...
         0];

% r_m2g = [r_m2*g*cos(r_q1+r_q2); -r_m2*g*sin(r_q1+r_q2); 0];
syms r_m2gi  % r_m2*g*cos(r_q1+r_q2)
syms r_m2gj  % -r_m2*g*sin(r_q1+r_q2)
r_m2g = [r_m2gi; r_m2gj; 0];

r_f2 = simplify(r_m2*r_h_ac2 - r_m2g);

r_dhc2 = [0; 0; r_I_G2z*(r_ddq1 + r_ddq2)]; % rate of change of angular momentum of body 2
r_tau2 = simplify(r_dhc2 - cross(r_f2, [r_lc2;0;0]) + tau_int);

%% 
r_h_ac1 = [-r_lc1*r_dq1^2; r_lc1*r_ddq1; 0];

% m1g = [m1*g*cos(q1); -m1*g*sin(q1); 0];
syms r_m1gi  % m1*g*cos(q1)
syms r_m1gj  % -m1*g*sin(q1)
r_m1g = [r_m1gi; r_m1gj; 0];

R12 = get_rotmax(R_T12);
r_f1 = r_m1*r_h_ac1 + R12*r_f2 - r_m1g;

r_dhc1 = [0; 0; r_I_G1z*r_ddq1]; % rate of change of angular momentum of body 1
r_tau1 = simplify(R12*r_tau2 - cross(r_f1, [r_lc1;0;0]) - cross(R12*r_f2, [r_l1-r_lc1;0;0]) + r_dhc1 + tau_int);

%%
% Human side
h_h_ac2 = [-h_lc2*(h_dq1 + h_dq2)^2 - h_l1*cos(h_q2)*h_dq1^2  + h_l1*sin(h_q2)*h_ddq2;...
         h_lc2*(h_ddq1 + h_ddq2) + h_l1*sin(h_q2)*h_dq1^2 + h_l1*cos(h_q2)*h_ddq1;...
         0];

% h_m2g = [h_m2*g*cos(h_q1+h_q2); -h_m2*g*sin(h_q1+h_q2); 0];
syms h_m2gi  % h_m2*g*cos(h_q1+h_q2)
syms h_m2gj  % -h_m2*g*sin(h_q1+h_q2)
h_m2g = [h_m2gi; h_m2gj; 0];

h_f2 = simplify(h_m2*h_h_ac2 - h_m2g);

h_dhc2 = [0; 0; h_I_G2z*(h_ddq1 + h_ddq2)]; % rate of change of angular momentum of body 2
h_tau2 = simplify(h_dhc2 - cross(h_f2, [h_lc2;0;0]) - tau_int);

%% 
h_h_ac1 = [-h_lc1*h_dq1^2; h_lc1*h_ddq1; 0];

% h_m1g = [h_m1*g*cos(h_q1); -h_m1*g*sin(h_q1); 0];
syms h_m1gi  % h_m1*g*cos(h_q1)
syms h_m1gj  % -h_m1*g*sin(h_q1)
h_m1g = [h_m1gi; h_m1gj; 0];

R12 = get_rotmax(H_T12);
h_f1 = h_m1*h_h_ac1 + R12*h_f2 - h_m1g;

h_dhc1 = [0; 0; h_I_G1z*h_ddq1]; % rate of change of angular momentum of body 1
h_tau1 = simplify(R12*h_tau2 - cross(h_f1, [h_lc1;0;0]) - cross(R12*h_f2, [h_l1-h_lc1;0;0]) + h_dhc1 - tau_int);


%% Collect terms
my_variables = [h_q1, h_dq1, h_ddq1, ...
                h_q2, h_dq2, h_ddq2, ...
                r_q1, r_dq1, r_ddq1, ...
                r_q2, r_dq2, r_ddq2];

r_tau1_eqn = subs(r_tau1, {h_q1, h_dq1, h_ddq1, r_q1, r_dq1, r_ddq1}, {pi/2, 0, 0, pi/2, 0, 0});
r_tau1_eqn = collect(r_tau1_eqn, my_variables)

r_tau2_eqn = subs(r_tau2, {h_q1, h_dq1, h_ddq1, r_q1, r_dq1, r_ddq1}, {pi/2, 0, 0, pi/2, 0, 0});
r_tau2_eqn = collect(r_tau2_eqn, my_variables)

h_tau1_eqn = subs(h_tau1, {h_q1, h_dq1, h_ddq1, r_q1, r_dq1, r_ddq1}, {pi/2, 0, 0, pi/2, 0, 0});
h_tau1_eqn = collect(h_tau1_eqn, my_variables)

h_tau2_eqn = subs(h_tau2, {h_q1, h_dq1, h_ddq1, r_q1, r_dq1, r_ddq1}, {pi/2, 0, 0, pi/2, 0, 0});
h_tau2_eqn = collect(h_tau2_eqn, my_variables)


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
