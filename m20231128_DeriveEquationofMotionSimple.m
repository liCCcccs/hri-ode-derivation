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

syms h_tau1 h_tau2 r_tau1 r_tau2

h_q1 = pi/2; h_dq1 = 0; h_ddq1 = 0;
r_q1 = pi/2; r_dq1 = 0; r_ddq1 = 0;
h_tau1 = 0; r_tau1 = 0;

% intermidiate variables
syms r_m2gi  % r_m2*g*cos(r_q1+r_q2)
syms r_m2gj  % -r_m2*g*sin(r_q1+r_q2)
syms r_m1gi  % m1*g*cos(q1)
syms r_m1gj  % -m1*g*sin(q1)
syms h_m2gi  % h_m2*g*cos(h_q1+h_q2)
syms h_m2gj  % -h_m2*g*sin(h_q1+h_q2)
syms h_m1gi  % h_m1*g*cos(h_q1)
syms h_m1gj  % -h_m1*g*sin(h_q1)
syms h_m1g h_m2g r_m1g r_m2g
    
eqns = [
    r_tau2 == (-Ks)*h_q2 + (-Ds)*h_dq2 + (As + Ks)*r_q2 + Ds*r_dq2 + (r_m2*r_lc2^2 + r_I_G2z)*r_ddq2 - r_lc2*r_m2gj

    h_tau2 == Ks*h_q2 + Ds*h_dq2 + (h_m2*h_lc2^2 + h_I_G2z)*h_ddq2 + (- As - Ks)*r_q2 + (-Ds)*r_dq2 - h_lc2*h_m2gj
    ];


%tobeElimitedVars = [];
%expr = eliminate(eqns, tobeElimitedVars);
ddthetaSolved = solve(eqns, [h_ddq2, r_ddq2]);
%ddthetaSolved = solve(expr, [ddq_2, ddq_3, ddq_4]);

simplify(ddthetaSolved.h_ddq2)
simplify(ddthetaSolved.r_ddq2)


vSym.h_m1 = 7.275;
vSym.h_m2 = 3.75;
vSym.r_m1 = 2;
vSym.r_m2 = 2;
vSym.g = 10;
vSym.h_I_G1z = 0.121; vSym.h_I_G2z = 0.055; vSym.r_I_G1z = 0.02; vSym.r_I_G2z = 0.02; 
vSym.h_l1 = 0.4; vSym.h_l2 = 0.4; vSym.r_l1 = 0.2; vSym.r_l2 = 0.2;
vSym.h_lc1 = 0.173; vSym.h_lc2 = 0.173; 
vSym.r_lc1 = 0.1; vSym.r_lc2 = 0.1; 
vSym.Ks = 10; vSym.Ds = 5; vSym.As = 0;
% vSym.K_AFz = 0; vSym.K_AFx = 0; vSym.K_AMy = 0;
% vSym.K_BFz = 0; vSym.K_BFx = 0; vSym.K_BMy = 0;
% vSym.D_AFz = 1000; vSym.D_AFx = -1000; vSym.D_AMy = 10;
% vSym.D_BFz = 0; vSym.D_BFx = 0; vSym.D_BMy = 0;
% vSym.K_AMy = 0; vSym.K_BMy = 0;
% vSym.D_AFz = 300; vSym.D_AFx = 300; vSym.D_AMy = 20;
% vSym.D_BFz = 300; vSym.D_BFx = 300; vSym.D_BMy = 20;

%% solve the dynamic equations
tspan = [0 1];
%tspan = 0:0.05:1;
y0 = [h_q1, h_dq1, -pi/8, 0, r_q1, r_dq1, -pi/8, 0];  %initial condition
tau = [0, 0, 0, 0];
[t,y] = ode23(@(t,y) runrobot(t,y, tau, vSym), tspan, y0);

%% Plot trajectory
% figure('Renderer', 'painters', 'Position', [300 300 800 800])
% plot(t, y(:, [3, 11]))
% legend(["h th", "r th"])


%% save data
total_steps = length(y);
hri_data = [y, zeros(total_steps, 4)];
header = ["h_q1", "h_dq1", "h_q2", "h_dq2", "r_q1", "r_dq1", "r_q2", "r_dq2", "tau1", "tau2", "tau3", "tau4"];
filename = 'hri_sim_data';
save_to_csv(hri_data, header, filename)
% save('hri_data.mat', 'hri_data')

%% Visualise
dt = 0.02;
h_q1 = y(:, 1); h_dq1 = y(:, 2);
h_q2 = y(:, 3); h_dq2 = y(:, 4);

r_q1 = y(:, 5); r_dq1 = y(:, 6);
r_q2 = y(:, 7); r_dq2 = y(:, 8);  

figure('Renderer', 'painters', 'Position', [300 300 600 500])
for i = 1:length(y)
    [h_j1, h_j2, h_j3, r_j1, r_j2, r_j3] = calc_joint_position(vSym, h_q1(i), h_q2(i), r_q1(i), r_q2(i));

    clf
    plot_3d_points({h_j1, h_j2, h_j3});
    plot_3d_points({r_j1, r_j2, r_j3});
    view(0, 90)
    legend("Human", "Robot")

    axis([-0.2,0.8,-0.6,0.4, -0.4,0.6])
    pause(dt)
end

%figure('Renderer', 'painters', 'Position', [300 300 800 700])
%plot(f_intB)

%% Calculate differetial equations
function dydt = runrobot(t, y, tau, vSym)

    vVar.h_q1 = y(1); vVar.h_dq1 = y(2); 
    vVar.h_q2 = y(3); vVar.h_dq2 = y(4); 
    vVar.r_q1 = y(5); vVar.r_dq1 = y(6); 
    vVar.r_q2 = y(7); vVar.r_dq2 = y(8); 

    vVar.h_tau1 = tau(1); vVar.r_tau1 = tau(3);
    if t < 4
        vVar.h_tau2 = 0;
        vVar.r_tau2 = 0 * sin(1.5*2*pi*t);
    else
        vVar.h_tau2 = 0;
        vVar.r_tau2 = 0;
    end

    %ddq1 = calc_equation_of_motion(ddthetaSolved.ddq1, symbols, vSym, vVar);
    h_ddq1 = 0; r_ddq1 = 0;
    [h_ddq2, r_ddq2] = calc_equation_of_motion(vSym, vVar);

    dydt = zeros(4,1);
    dydt(1) = vVar.h_dq1; dydt(2) = h_ddq1;
    dydt(3) = vVar.h_dq2; dydt(4) = h_ddq2;
    dydt(5) = vVar.r_dq1; dydt(6) = r_ddq1;
    dydt(7) = vVar.r_dq2; dydt(8) = r_ddq2;
end

function [h_ddq2, r_ddq2] = calc_equation_of_motion(vSym, vVar)
   

    % Assign values to symbols
    h_m1 = vSym.h_m1; h_m2 = vSym.h_m2;
    r_m1 = vSym.r_m1; r_m2 = vSym.r_m2;
    g = vSym.g;
    h_I_G1z = vSym.h_I_G1z; h_I_G2z = vSym.h_I_G2z; 
    r_I_G1z = vSym.r_I_G1z; r_I_G2z = vSym.r_I_G2z; 

    h_l1 = vSym.h_l1; h_l2 = vSym.h_l2; 
    r_l1 = vSym.r_l1; r_l2 = vSym.r_l2;
    h_lc1 = vSym.h_lc1; h_lc2 = vSym.h_lc2; 
    r_lc1 = vSym.r_lc1; r_lc2=  vSym.r_lc2; 
    Ks = vSym.Ks; Ds = vSym.Ds; As = vSym.As;

    h_q1 = vVar.h_q1; h_dq1 = vVar.h_dq1;
    h_q2 = vVar.h_q2;  h_dq2 = vVar.h_dq2;
    r_q1 = vVar.r_q1;  r_dq1 = vVar.r_dq1;
    r_q2 = vVar.r_q2;  r_dq2 = vVar.r_dq2;

    h_tau1 = vVar.h_tau1; h_tau2 = vVar.h_tau2; 
    r_tau1 = vVar.r_tau1; r_tau2 = vVar.r_tau2; 

    r_m2gi = r_m2*g*cos(r_q1+r_q2);
    r_m2gj = -r_m2*g*sin(r_q1+r_q2);
    r_m1gi = r_m1*g*cos(r_q1);
    r_m1gj = -r_m1*g*sin(r_q1);
    h_m2gi = h_m2*g*cos(h_q1+h_q2);
    h_m2gj = -h_m2*g*sin(h_q1+h_q2);
    h_m1gi = h_m1*g*cos(h_q1);
    h_m1gj = -h_m1*g*sin(h_q1);

    h_ddq2 = (h_tau2 - Ds*h_dq2 + As*r_q2 - Ks*h_q2 + Ds*r_dq2 + Ks*r_q2 + h_lc2*h_m2gj)/(h_m2*h_lc2^2 + h_I_G2z);
    r_ddq2 = (r_tau2 + Ds*h_dq2 - As*r_q2 + Ks*h_q2 - Ds*r_dq2 - Ks*r_q2 + r_lc2*r_m2gj)/(r_m2*r_lc2^2 + r_I_G2z);

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


function [h_j1, h_j2, h_j3, r_j1, r_j2, r_j3] = calc_joint_position(vSym, h_q1, h_q2, r_q1, r_q2)
    H_T01 = calc_homotransmtx(0, 0, 0, -pi/2 + h_q1);
    H_T12 = calc_homotransmtx(vSym.h_l1, 0, 0, h_q2);
    H_T23 = calc_homotransmtx(vSym.h_l2, 0, 0, 0);
    H_T02 = H_T01 * H_T12;
    H_T03 = H_T02 * H_T23;

    R_T01 = calc_homotransmtx(0, 0, 0, -pi/2 + r_q1);
    R_T12 = calc_homotransmtx(vSym.r_l1, 0, 0, r_q2);
    R_T23 = calc_homotransmtx(vSym.r_l2, 0, 0, 0);
    R_T02 = R_T01 * R_T12;
    R_T03 = R_T02 * R_T23;

    h_j1 = [0; 0; 0];
    h_j2 = H_T02(1:3, 4);
    h_j3 = H_T03(1:3, 4);
    
    r_j1 = [0; 0; 0] + [vSym.h_l1 - vSym.r_l1; 0; 0];
    r_j2 = R_T02(1:3, 4) + [vSym.h_l1 - vSym.r_l1; 0; 0];
    r_j3 = R_T03(1:3, 4) + [vSym.h_l1 - vSym.r_l1; 0; 0];
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
    %scatter3(xs, ys, zs, 'o', 'LineWidth', 2); % 'o' is for circle marker

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

function save_to_csv(data, header, filename)
    T = array2table(data);
    T.Properties.VariableNames(1:length(header)) = header;
    writetable(T,[filename, '.csv']);
    %print(['saved to ', filename]);
end
