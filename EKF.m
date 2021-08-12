function [n_est, d_est, alphaij_est, P_new] = EKF(x_est_prev,U,P_prev,z_measured) %inputs at time t-1

global G_fit_det outflow1STD outflow2STD outflow3STD outflow4STD std_acc_meas std_transferflow_meas std_randomwalk_demand std_randomwalk_alpha

Q = diag([outflow1STD^2 outflow2STD^2 outflow3STD^2 outflow4STD^2 std_randomwalk_demand^2*ones(1,4) std_randomwalk_alpha^2*ones(1,6)]);%E{ksi*ksi} ---> ksi = [error in MFDs ; error in random walks]'---> size=(14*1);
R = [Q(1:4,1:4) zeros(4,10);zeros(4,4) (std_acc_meas^2)*eye(4) zeros(4,6);zeros(6,8) (std_transferflow_meas^2)*eye(6)];%E{eta*eta}
% M = [Q(1:4,1:4) zeros(4,10);zeros(10,4) zeros(10,10)];%E{ksi*eta} -> Related to Markos definition
M = [zeros(4,4) zeros(4,10);zeros(10,4) zeros(10,10)];%E{ksi*eta} -> Related to Markos definition


n_est_prev = x_est_prev(1:4);
d_est_prev = x_est_prev(5:8);
alphaij_est_prev = x_est_prev(9:14);
U14 = U(1); U24 = U(2); U34 = U(3); U41 = U(4); U42 = U(5); U43 = U(6);

[n_est_prior , d_est_prior, alphaij_est_prior]= plant4_modified_new(d_est_prev,alphaij_est_prev,n_est_prev,[U14;U24;U34;U41;U42;U43]);
x_est_prior = [n_est_prior;d_est_prior;alphaij_est_prior];
A = lin_form_new(x_est_prev, [U14;U24;U34;U41;U42;U43]);%df/dx
C = measure_gain_new(x_est_prev,[U14;U24;U34;U41;U42;U43]);%dg/dx
W = matrix_gama_new(x_est_prev,[U14;U24;U34;U41;U42;U43]);%df/dksi
S = matrix_sigma_new(x_est_prev,[U14;U24;U34;U41;U42;U43]);%dg/deta
P_pred = A*P_prev*A'+W*Q*W';
K = (P_pred*C'+W*M*S')*(C*P_pred*C'+S*R*S'+S*M'*W'*C'+C*W'*M*S')^-1;
P_new = P_pred - K*(C*P_pred-S*M'*W');

%Correction step in Kalman filter
x_est_new = x_est_prior + K * (z_measured - [x_est_prior(1);x_est_prior(2);x_est_prior(3);x_est_prior(4);...
                                            U14*(1-x_est_prior(9))*polyval(G_fit_det(1,:),x_est_prior(1));...
                                            U24*(1-x_est_prior(10))*polyval(G_fit_det(2,:),x_est_prior(2));...
                                            U34*(1-x_est_prior(11))*polyval(G_fit_det(3,:),x_est_prior(3));...
                                            U41*x_est_prior(12)*polyval(G_fit_det(4,:),x_est_prior(4));...
                                            U42*x_est_prior(13)*polyval(G_fit_det(4,:),x_est_prior(4));...
                                            U43*x_est_prior(14)*polyval(G_fit_det(4,:),x_est_prior(4));]);                                       
n_est = x_est_new(1:4);
d_est = x_est_new(5:8);
alphaij_est = x_est_new(9:14);
