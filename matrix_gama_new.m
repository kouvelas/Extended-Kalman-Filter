function G = matrix_gama_new(x_est,U)

global G_fit_det

n1 = x_est(1);  n2 = x_est(2);  n3 = x_est(3);  n4 = x_est(4);

U14 = U(1); U24 = U(2); U34 = U(3); U41 = U(4); U42 = U(5); U43 = U(6);

a11 = x_est(9);  a22 = x_est(10);  a33 = x_est(11);  a41 = x_est(12);  a42 = x_est(13);  a43 = x_est(14); 

G = zeros(14,14);

%n1(k+1) = n1(k) + d1(k) - a11(k)*G1(n1) + U41(k)*a41*G4(n4) - (1-a11(k))*G1(n1)*U14;
G(1,1) = -a11*polyval(G_fit_det(1,:),n1) - (1-a11)*polyval(G_fit_det(1,:),n1)*U14;
G(1,4) = U41*a41*polyval(G_fit_det(4,:),n4);

%n2(k+1) = n2(k) + d2(k) - a22(k)*G2(n2) + U42(k)*a42*G4(n4) - (1-a22(k))*G2(n2)*U24;
G(2,2) = -a22*polyval(G_fit_det(2,:),n2)-(1-a22)*polyval(G_fit_det(2,:),n2)*U24;
G(2,4) = U42*a42*polyval(G_fit_det(4,:),n4);

%n3(k+1) = n3(k) + d3(k) - a33(k)*G3(n3) + U43(k)*a43*G4(n4) - (1-a33(k))*G3(n3)*U34;
G(3,3) = -a33*polyval(G_fit_det(3,:),n3)-(1-a33)*polyval(G_fit_det(3,:),n3)*U34;
G(3,4) = U43*a43*polyval(G_fit_det(4,:),n4);

%n4(k+1) = n4(k) + d4(k) - (1-a41(k)-a42(k)-a43(k))*G4(n4) - a41*U41*G4(n4) - a42*U42*G4(n4) - a43*U43*G4(n4) + (1-a11)*U14*G1(n1) + (1-a22)*U24*G2(n2) + (1-a33)*U34*G3(n3);
G(4,1) = (1-a11)*U14*polyval(G_fit_det(1,:),n1);
G(4,2) = (1-a22)*U24*polyval(G_fit_det(2,:),n2);
G(4,3) = (1-a33)*U34*polyval(G_fit_det(3,:),n3);
G(4,4) = -(1-a41-a42-a43)*polyval(G_fit_det(4,:),n4) - a41*U41*polyval(G_fit_det(4,:),n4) - a42*U42*polyval(G_fit_det(4,:),n4) - a43*U43*polyval(G_fit_det(4,:),n4);


G(5:14,:)=[zeros(10,4) eye(10)];
