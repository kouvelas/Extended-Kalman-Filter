clear all
close all
clc
load G_fit_det_1
load G_fit_det_2
load P_fit_90
a1 = 1.15;
P_fit_90(1,:) = [P_fit_90(1,1)*a1^3 P_fit_90(1,2)*a1^2 P_fit_90(1,3)*a1 P_fit_90(1,4)];
P_fit = P_fit_90;
load L_ave
global G_fit_det std_acc_meas std_transferflow_meas outflow1STD outflow2STD outflow3STD outflow4STD n_max std_randomwalk_demand  std_randomwalk_alpha
G_fit_det = G_fit_det_1;
outflow1STD = [0.01, 0.02, 0.03];
outflow2STD = [0.01, 0.02, 0.03];
outflow3STD = [0.01, 0.02, 0.03];
outflow4STD = [0.01, 0.02, 0.03];
%Measurement noise
std_acc_meas = [0.01, 0.02, 0.03];
std_transferflow_meas = [2, 3, 4];
std_randomwalk_demand = [30, 40, 50];
std_randomwalk_alpha = [0.01, 0.02, 0.03];
n_max = [15390,6210,6480,19150];

all_comb = {outflow1STD, outflow2STD, outflow3STD, outflow4STD, std_acc_meas, std_transferflow_meas, std_randomwalk_demand, std_randomwalk_alpha};
NN = length(all_comb);
[all_comb{:}] = ndgrid(all_comb{:});
all_comb = reshape(cat(NN+1,all_comb{:}),[],NN);
size(all_comb)

tasos_res=[];

for my_row=1:6561
    
tic

load G_fit_det_1
load G_fit_det_2
load P_fit_90
a1 = 1.15;
P_fit_90(1,:) = [P_fit_90(1,1)*a1^3 P_fit_90(1,2)*a1^2 P_fit_90(1,3)*a1 P_fit_90(1,4)];
P_fit = P_fit_90;
load L_ave
global G_fit_det std_acc_meas std_transferflow_meas outflow1STD outflow2STD outflow3STD outflow4STD n_max std_randomwalk_demand  std_randomwalk_alpha
G_fit_det = G_fit_det_1;
n_max = [15390,6210,6480,19150];

outflow1STD = all_comb(my_row,1);
outflow2STD = all_comb(my_row,2);
outflow3STD = all_comb(my_row,3);
outflow4STD = all_comb(my_row,4);
%Measurement noise
std_acc_meas = all_comb(my_row,5);
std_transferflow_meas = all_comb(my_row,6);
std_randomwalk_demand = all_comb(my_row,7);
std_randomwalk_alpha = all_comb(my_row,8);


%%
% Plotting the MFDs
load('aimDatoldaimsun.txt');
aimDat=aimDatoldaimsun;
nu1 = std_acc_meas * randn(size(aimDat,1),4);
nu2 = std_transferflow_meas * randn(size(aimDat,1),6);

%%
r_init = 3;% the first time interval to estimate
t_end = 180;%size(aimDat,1);

for r=1:t_end
    if (r <= 82)
        G_fit_det = G_fit_det_1;
    else
        G_fit_det = G_fit_det_2;
    end
    if r < r_init
        P_prev = eye(14);
        save Pmatrix.txt P_prev -ascii;
    else
        P_prev = load('Pmatrix.txt');
        load('ekfDat_test.txt');
        x_prev = ekfDat_test(r-1,1:14);
        z_measured = [(1+nu1(r,:)).* aimDat(r,19:22) aimDat(r+1,[10 12 14 15 16 17])+nu2(r+1,:)]';
        U = ones(6,1);
        [n_est_new, d_est_new, alphaij_est_new, P_new] = EKF(x_prev,U,P_prev,z_measured);
        save Pmatrix.txt P_new -ascii;
    end
    
    if r < r_init
        if (r==1)
            dat = [aimDat(r,19:22) aimDat(r+1,5:8) aimDat(r+1,9)/(aimDat(r+1,9)+aimDat(r+1,10)) aimDat(r+1,11)/(aimDat(r+1,11)+aimDat(r+1,12))...
                                               aimDat(r+1,13)/(aimDat(r+1,13)+aimDat(r+1,14)) aimDat(r+1,15)/(sum(aimDat(r+1,15:18))) ...
                                               aimDat(r+1,16)/(sum(aimDat(r+1,15:18))) aimDat(r+1,17)/(sum(aimDat(r+1,15:18)))];
            save ekfDat_test.txt dat -ascii;
        else
            load('ekfDat_test.txt');
            dat = [ekfDat_test;aimDat(r,19:22) aimDat(r+1,5:8) aimDat(r+1,9)/(aimDat(r+1,9)+aimDat(r+1,10)) aimDat(r+1,11)/(aimDat(r+1,11)+aimDat(r+1,12))...
                                               aimDat(r+1,13)/(aimDat(r+1,13)+aimDat(r+1,14)) aimDat(r+1,15)/(sum(aimDat(r+1,15:18))) ...
                                               aimDat(r+1,16)/(sum(aimDat(r+1,15:18))) aimDat(r+1,17)/(sum(aimDat(r+1,15:18)))];
            save ekfDat_test.txt dat -ascii;
        end
    else
        dat = [ekfDat_test; n_est_new' d_est_new' alphaij_est_new'];
        save ekfDat_test.txt dat -ascii;
    end
end

clearvars -except r_init t_end aimDat nu1 my_row tasos_res all_comb
close all
clc
%%
load ekfDat_test.txt
ekfDat = ekfDat_test;

X1 = [aimDat(r_init-1:t_end,19)+nu1(r_init-1:t_end,1),aimDat(r_init-1:t_end,20)+nu1(r_init-1:t_end,2),aimDat(r_init-1:t_end,21)+nu1(r_init-1:t_end,3),aimDat(r_init-1:t_end,22)+nu1(r_init-1:t_end,4)];
X2 = [ekfDat(r_init-1:t_end,1),ekfDat(r_init-1:t_end,2),ekfDat(r_init-1:t_end,3),ekfDat(r_init-1:t_end,4)];
PI_n1 = sqrt(1/(179*4)*sum((X1(:,1)-X2(:,1)).^2));
PI2_n1 = PI_n1/(sqrt(1/(179*4))*sum(X1(:,1)));
PI_n2 = sqrt(1/(179*4)*sum((X1(:,2)-X2(:,2)).^2));
PI2_n2 = PI_n2/(sqrt(1/(179*4))*sum(X1(:,2)));
PI_n3 = sqrt(1/(179*4)*sum((X1(:,3)-X2(:,3)).^2));
PI2_n3 = PI_n3/(sqrt(1/(179*4))*sum(X1(:,3)));
PI_n4 = sqrt(1/(179*4)*sum((X1(:,4)-X2(:,4)).^2));
PI2_n4 = PI_n4/(sqrt(1/(179*4))*sum(X1(:,4)));

X1 = [aimDat(r_init-1:t_end,5),aimDat(r_init-1:t_end,6),aimDat(r_init-1:t_end,7),aimDat(r_init-1:t_end,8)];
X2 = [ekfDat(r_init-1:t_end,5),ekfDat(r_init-1:t_end,6),ekfDat(r_init-1:t_end,7),ekfDat(r_init-1:t_end,8)];
PI_d1 = sqrt(1/(179*4)*sum((X1(:,1)-X2(:,1)).^2));
PI2_d1 = PI_d1/(sqrt(1/(179*4))*sum(X1(:,1)));
PI_d2 = sqrt(1/(179*4)*sum((X1(:,2)-X2(:,2)).^2));
PI2_d2 = PI_d2/(sqrt(1/(179*4))*sum(X1(:,2)));
PI_d3 = sqrt(1/(179*4)*sum((X1(:,3)-X2(:,3)).^2));
PI2_d3 = PI_d3/(sqrt(1/(179*4))*sum(X1(:,3)));
PI_d4 = sqrt(1/(179*4)*sum((X1(:,4)-X2(:,4)).^2));
PI2_d4 = PI_d4/(sqrt(1/(179*4))*sum(X1(:,4)));

xx1 = aimDat(r_init-1:t_end,9)./(aimDat(r_init-1:t_end,9)+aimDat(r_init-1:t_end,10));
xx1(isnan(xx1)) = 1;
xx2 = aimDat(r_init-1:t_end,11)./(aimDat(r_init-1:t_end,11)+aimDat(r_init-1:t_end,12));
xx2(isnan(xx2)) = 1;
xx3 = aimDat(r_init-1:t_end,13)./(aimDat(r_init-1:t_end,13)+aimDat(r_init-1:t_end,14));
xx3(isnan(xx3)) = 1;
xx4 = aimDat(r_init-1:t_end,15)./sum(aimDat(r_init-1:t_end,15:18),2);
xx4(isnan(xx4)) = 1;
xx5 = aimDat(r_init-1:t_end,16)./sum(aimDat(r_init-1:t_end,15:18),2);
xx5(isnan(xx5)) = 1;
xx6 = aimDat(r_init-1:t_end,17)./sum(aimDat(r_init-1:t_end,15:18),2);
xx6(isnan(xx6)) = 1;
X1 = [xx1,xx2,xx3,xx4,xx5,xx6];
X2 = [ekfDat(r_init-1:t_end,9),ekfDat(r_init-1:t_end,10),ekfDat(r_init-1:t_end,11),ekfDat(r_init-1:t_end,12),ekfDat(r_init-1:t_end,13),ekfDat(r_init-1:t_end,14)];
PI_a1 = sqrt(1/(179*6)*sum((X1(:,1)-X2(:,1)).^2));
PI2_a1 = PI_a1/(sqrt(1/(179*6))*sum(X1(:,1)));
PI_a2 = sqrt(1/(179*6)*sum((X1(:,2)-X2(:,2)).^2));
PI2_a2 = PI_a2/(sqrt(1/(179*6))*sum(X1(:,2)));
PI_a3 = sqrt(1/(179*6)*sum((X1(:,3)-X2(:,3)).^2));
PI2_a3 = PI_a3/(sqrt(1/(179*6))*sum(X1(:,3)));
PI_a4 = sqrt(1/(179*6)*sum((X1(:,4)-X2(:,4)).^2));
PI2_a4 = PI_a4/(sqrt(1/(179*6))*sum(X1(:,4)));
PI_a5 = sqrt(1/(179*6)*sum((X1(:,5)-X2(:,5)).^2));
PI2_a5 = PI_a5/(sqrt(1/(179*6))*sum(X1(:,5)));
PI_a6 = sqrt(1/(179*6)*sum((X1(:,6)-X2(:,6)).^2));
PI2_a6 = PI_a6/(sqrt(1/(179*6))*sum(X1(:,6)));

new_row = [my_row PI_n1 PI2_n1 PI_n2 PI2_n2 PI_n3 PI2_n3 PI_n4 PI2_n4];
new_row = [new_row PI_d1 PI2_d1 PI_d2 PI2_d2 PI_d3 PI2_d3 PI_d4 PI2_d4];
new_row = [new_row PI_a1 PI2_a1 PI_a2 PI2_a2 PI_a3 PI2_a3 PI_a4 PI2_a4 PI_a5 PI2_a5 PI_a6 PI2_a6]

tasos_res=[tasos_res; new_row];

end

save tasos_res2.txt tasos_res -ascii;

%%
clearvars -except r_init t_end aimDat nu1
close all
clc
%%
load ekfDat_test.txt
ekfDat = ekfDat_test;

toc

%%
%plotting accumulation (real vs estimated)
figure

subplot(2,2,1)
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,19)+nu1(r_init-1:t_end,1),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,1),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('n_1')
legend('real','est.','Location','Best')


subplot(2,2,2)
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,20)+nu1(r_init-1:t_end,2),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,2),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('n_2')
legend('real','est.','Location','Best')

subplot(2,2,3)
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,21)+nu1(r_init-1:t_end,3),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,3),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('n_3')
legend('real','est.','Location','Best')

subplot(2,2,4)
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,22)+nu1(r_init-1:t_end,4),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,4),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('n_4')
legend('real','est.','Location','Best')

%%
%plotting demands (real vs estimated)
figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,5),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,5),'r','Linewidth',1.5)
grid
ylim([0 1300])
box on
xlabel('time')
ylabel('d_1')
legend('real','est.','Location','Best')


figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,6),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,6),'r','Linewidth',1.5)
grid
ylim([0 1300])
box on
xlabel('time')
ylabel('d_2')
legend('real','est.','Location','Best')

figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,7),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,7),'r','Linewidth',1.5)
grid
ylim([0 1300])
box on
xlabel('time')
ylabel('d_3')
legend('real','est.','Location','Best')

figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,8),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,8),'r','Linewidth',1.5)
grid
ylim([0 1300])
box on
xlabel('time')
ylabel('d_4')
legend('real','est.','Location','Best')


%%
%plotting alphas (real vs estimated)
figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,9)./(aimDat(r_init-1:t_end,9)+aimDat(r_init-1:t_end,10)),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,9),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('\alpha_{11}')
legend('real','est.','Location','Best')
ylim([0 1])


figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,11)./(aimDat(r_init-1:t_end,11)+aimDat(r_init-1:t_end,12)),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,10),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('\alpha_{22}')
legend('real','est.','Location','Best')
ylim([0 1])


figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,13)./(aimDat(r_init-1:t_end,13)+aimDat(r_init-1:t_end,14)),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,11),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('\alpha_{33}')
legend('real','est.','Location','Best')
ylim([0 1])

figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,15)./sum(aimDat(r_init-1:t_end,15:18),2),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,12),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('\alpha_{41}')
legend('real','est.','Location','Best')
ylim([0 1])

figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,16)./sum(aimDat(r_init-1:t_end,15:18),2),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,13),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('\alpha_{42}')
legend('real','est.','Location','Best')
ylim([0 1])


figure
hold on
plot(r_init-1:t_end,aimDat(r_init-1:t_end,17)./sum(aimDat(r_init-1:t_end,15:18),2),'b-.','Linewidth',1.5)
plot(r_init-1:t_end,ekfDat(r_init-1:t_end,14),'r','Linewidth',1.5)
grid
xlabel('time')
box on
ylabel('\alpha_{43}')
legend('real','est.','Location','Best')
ylim([0 1])

