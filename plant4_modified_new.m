function [n_next , d_next , alphaij_next]= plant4_modified_new(d,alpha,n_current,u)

global G_fit_det n_max

N1_current=n_current(1);N2_current=n_current(2);N3_current=n_current(3);N4_current=n_current(4);

U14=u(1);U24=u(2);U34=u(3);
U41=u(4);U42=u(5);U43=u(6);

d1=d(1);d2=d(2);d3=d(3);d4=d(4);

a11 = alpha(1); a22 = alpha(2); a33 = alpha(3);
a41 = alpha(4); a42 = alpha(5); a43 = alpha(6);   
%%
%updating equation
%Region #1
if (N1_current <= n_max(1))
    M11 = polyval(G_fit_det(1,:),N1_current)*a11;
    M1R = polyval(G_fit_det(1,:),N1_current)*(1-a11);
else
    M11 = 0;
    M1R = 0;
end
%Region #2
if (N2_current <= n_max(2))
    M22 = polyval(G_fit_det(2,:),N2_current)*a22;
    M2R = polyval(G_fit_det(2,:),N2_current)*(1-a22);
else
    M2R = 0;
    M22 = 0;
end
%Region #3
if (N3_current <= n_max(3))
    M33 = polyval(G_fit_det(3,:),N3_current)*a33;
    M3R = polyval(G_fit_det(3,:),N3_current)*(1-a33);
else
    M33 = 0;
    M3R = 0;
end
%Region #4
if (N4_current <= n_max(4))
    M41 = polyval(G_fit_det(4,:),N4_current)*a41;
    M42 = polyval(G_fit_det(4,:),N4_current)*a42;
    M43 = polyval(G_fit_det(4,:),N4_current)*a43;
    M44 = polyval(G_fit_det(4,:),N4_current)*(1-a41-a42-a43);
else
    M41 = 0;
    M42 = 0;
    M43 = 0;
    M44 = 0;
end
%Region #1
N1_next = N1_current + d1 - M11 + U41*M41 - U14*M1R;
%Region #2
N2_next = N2_current + d2 - M22 + U42*M42 - U24*M2R;
%Region #3
N3_next = N3_current + d3 - M33 + U43*M43 - U34*M3R;
%Region #4
N4_next = N4_current + d4 - M44 - U41*M41 - U42*M42 - U43*M43 + U14*M1R + U24*M2R + U34*M3R;


n_next=zeros(4,1);
n_next(1) = max(0,N1_next) ; n_next(2) = max(0,N2_next) ; n_next(3) = max(0,N3_next) ; n_next(4) = max(0,N4_next);

d_next=zeros(4,1);
d_next(1,1) = d1; d_next(2,1) = d2; d_next(3,1) = d3; d_next(4,1) = d4;

alphaij_next=zeros(6,1);
alphaij_next(1,1) = a11; alphaij_next(2,1) = a22; alphaij_next(3,1) = a33; 
alphaij_next(4,1) = a41; alphaij_next(5,1) = a42; alphaij_next(6,1) = a43; 

