function C = measure_gain_new(x_est,U)
global G_fit_det

n1=x_est(1);n2=x_est(2);n3=x_est(3);n4=x_est(4);

U14 = U(1); U24 = U(2); U34 = U(3); U41 = U(4); U42 = U(5); U43 = U(6);

a11=x_est(9);a22=x_est(10);a33=x_est(11);a41=x_est(12);a42=x_est(13);a43=x_est(14);

C = zeros(6,14);


% M1R = (1-a11)*G1(n1)*U14;
C(1,1) = (1-a11)*U14*polyval([3*G_fit_det(1,1) 2*G_fit_det(1,2) G_fit_det(1,3)],n1);
C(1,9) = -U14*polyval(G_fit_det(1,:),n1);

% M2R = (1-a22)*G2(n2)*U24;
C(2,2) = (1-a22)*U24*polyval([3*G_fit_det(2,1) 2*G_fit_det(2,2) G_fit_det(2,3)],n2);
C(2,10) = -U24*polyval(G_fit_det(2,:),n2);


% M3R = (1-a33)*G3(n3)*U34;
C(3,3) = (1-a33)*U34*polyval([3*G_fit_det(3,1) 2*G_fit_det(3,2) G_fit_det(3,3)],n3);
C(3,11) = -U34*polyval(G_fit_det(3,:),n3);


% M41 = a41*G4(n4)*U41;
C(4,4) = a41*U41*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);
C(4,12) = U41*polyval(G_fit_det(4,:),n4);

% M42 = a42*G4(n4)*U42;
C(5,4) = a42*U42*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);
C(5,13) = U42*polyval(G_fit_det(4,:),n4);

% M43 = a43*G4(n4)*U43;
C(6,4) = a43*U43*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);
C(6,14) = U43*polyval(G_fit_det(4,:),n4);


temp = zeros(10,14);
temp(1:4,:) = [1 zeros(1,13);0 1 zeros(1,12); 0 0 1 zeros(1,11); 0 0 0 1 zeros(1,10)];
temp(5:10,:) = C(1:6,:);
C = temp;
