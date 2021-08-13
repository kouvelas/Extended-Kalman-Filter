function A = lin_form_new(x_est,U)      
%state variables=[n1 n2 n3 n4 d1 d2 d3 d4 a11 a22 a33 a41 a42 a43]
global G_fit_det

n1=x_est(1) ; n2=x_est(2) ; n3=x_est(3) ; n4=x_est(4);
a11=x_est(9) ; a22=x_est(10) ; a33=x_est(11) ; a41=x_est(12) ; a42=x_est(13) ; a43=x_est(14);     
U14 = U(1); U24 = U(2); U34 = U(3); U41 = U(4); U42 = U(5); U43 = U(6);

A = zeros(14,14);

%%
%Region 1
%n1(k+1) = n1(k) + d1(k) - a11(k)*G1(n1) + U41(k)*a41*G4(n4) - U14(k)*(1-a11)*G1(n1);
A(1,1) = 1 - (a11+U14*(1-a11))*polyval([3*G_fit_det(1,1) 2*G_fit_det(1,2) G_fit_det(1,3)],n1);%w.r.t n1
A(1,4) = a41*U41*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);%w.r.t n4
A(1,5) = 1;%w.r.t d1
A(1,9) = -polyval(G_fit_det(1,:),n1)+U14*polyval(G_fit_det(1,:),n1);%w.r.t a11
A(1,12) = U41*polyval(G_fit_det(4,:),n4);%w.r.t a41

%Region 2
%n2(k+1) = n2(k) + d2(k) - a22(k)*G2(n2) + U42(k)*a42*G4(n4) - U24(k)*(1-a22)*G2(n2);
A(2,2) = 1 - (a22+U24*(1-a22))*polyval([3*G_fit_det(2,1) 2*G_fit_det(2,2) G_fit_det(2,3)],n2);%w.r.t n2
A(2,4) = a42*U42*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);%w.r.t n4
A(2,6) = 1;%w.r.t d2
A(2,10) = -polyval(G_fit_det(2,:),n2)+U24*polyval(G_fit_det(2,:),n2);%w.r.t a22
A(2,13) = U42*polyval(G_fit_det(4,:),n4);%w.r.t a42


%Region 3
%n3(k+1) = n3(k) + d3(k) - a33(k)*G3(n3) + U43(k)*a43*G4(n4) - U34(k)*(1-a33)*G2(n3);
A(3,3) = 1 - (a33+U34*(1-a33))*polyval([3*G_fit_det(3,1) 2*G_fit_det(3,2) G_fit_det(3,3)],n3);%w.r.t n3
A(3,4) = a43*U43*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);%w.r.t n4
A(3,7) = 1;%w.r.t d3
A(3,11) = -polyval(G_fit_det(3,:),n3)+U34*polyval(G_fit_det(3,:),n3);%w.r.t a33
A(3,14) = U43*polyval(G_fit_det(4,:),n4);%w.r.t a43


%Region 4
%n4(k+1) = n4(k) + d4(k) - a41(k)*G4(n4)*U41 - a42(k)*G4(n4)*U42 - a43(k)*G4(n4)*U43 - (1-a41-a42-a43)*G4(n4) + (1-a11(k))*G1(n1)*U14 + (1-a22(k))*G2(n2)*U24 + (1-a33(k))*G3(n3)*U34;
A(4,1)  = (1-a11)*U14*polyval([3*G_fit_det(1,1) 2*G_fit_det(1,2) G_fit_det(1,3)],n1);%w.r.t n1
A(4,2)  = (1-a22)*U24*polyval([3*G_fit_det(2,1) 2*G_fit_det(2,2) G_fit_det(2,3)],n2);%w.r.t n2
A(4,3)  = (1-a33)*U34*polyval([3*G_fit_det(3,1) 2*G_fit_det(3,2) G_fit_det(3,3)],n3);%w.r.t n3
A(4,4)  = 1 - (a41*U41+a42*U42+a43*U43+1-a41-a42-a43)*polyval([3*G_fit_det(4,1) 2*G_fit_det(4,2) G_fit_det(4,3)],n4);%w.r.t n4
A(4,8)  = 1;%w.r.t d14
A(4,9)  = -U14*polyval(G_fit_det(1,:),n1);%w.r.t a11
A(4,10) = -U24*polyval(G_fit_det(2,:),n2);%w.r.t a22
A(4,11) = -U34*polyval(G_fit_det(3,:),n3);%w.r.t a33
A(4,12) = -U41*polyval(G_fit_det(4,:),n4)+polyval(G_fit_det(4,:),n4);%w.r.t a41
A(4,13) = -U42*polyval(G_fit_det(4,:),n4)+polyval(G_fit_det(4,:),n4);%w.r.t a42
A(4,14) = -U43*polyval(G_fit_det(4,:),n4)+polyval(G_fit_det(4,:),n4);%w.r.t a43


for i=5:14
   A(i,i) = 1; 
end



