function S = matrix_sigma_new(x_est,U)


global G_fit_det

n1=x_est(2);n2=x_est(2);n3=x_est(3);n4=x_est(4);

U14 = U(1); U24 = U(2); U34 = U(3); U41 = U(4); U42 = U(5); U43 = U(6);

a11=x_est(9);a22=x_est(10);a33=x_est(11);a41=x_est(12);a42=x_est(13);a43=x_est(14);

S = zeros(6,14);

% M1R = (1+ksi1) * (1-a11) * G1(n1) * U14 + gama5;
S(1,1) = (1-a11)*U14*polyval(G_fit_det(1,:),n1); % dM14/dksi1
S(1,9) = 1;


% M2R = (1+ksi2) * (1-a22) * G2(n2) * U24 + gama6;
S(2,2) = (1-a22)*U24*polyval(G_fit_det(2,:),n2); % dM24/dksi2
S(2,10) = 1;


% M3R = (1+ksi2) * (1-a33) * G3(n3) * U34 + gama7;
S(3,3) = (1-a33)*U34*polyval(G_fit_det(3,:),n3);
S(3,11) = 1;


% M41 = (1+ksi4) * a41 * G4(n4) * U41 + gama8;
S(4,4) = a41*U41*polyval(G_fit_det(4,:),n4);
S(4,12) = 1;


% M42 = (1+ksi4) * a42 * G4(n4) * U42 + gama9;
S(5,4) = a42*U42*polyval(G_fit_det(4,:),n4);
S(5,13) = 1;


% M43 = (1+ksi4) * a43 * G4(n4) * U43 + gama10;
S(6,4) = a43*U43*polyval(G_fit_det(4,:),n4);
S(6,14) = 1;


temp = zeros(10,14);
temp(1:4,:) = [zeros(1,4) n1 zeros(1,9); zeros(1,5) n2 zeros(1,8); zeros(1,6) n3 zeros(1,7); zeros(1,7) n4 zeros(1,6);];
temp(5:10,:) = S(1:6,:);
S = temp;
