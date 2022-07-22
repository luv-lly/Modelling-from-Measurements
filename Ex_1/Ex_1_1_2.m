%% opt-DMD
clear all; close all; clc

%Datasets

load HL_population.mat

train = 0.8*length(HL);
HL_train = HL(:,1:train);
t_train = time(1:train);

[w,e,b] = optdmd(HL_train,t_train,2,1);

X = w*diag(b)*exp(e*time);     

figure(1)

h = plot (years,H,'r--o',years,L,'b--o',years, abs(X(1,:)),'r',years, abs(X(2,:)),'b');
hold on
fill([years(end)-6 years(end) years(end) years(end)-6],[0 0 max(H) max(H)],'c','FaceAlpha',0.1,'LineStyle','-')
grid on
xlim([1845 1903])
ylabel ('Population State');
xlabel ('Years');
title ('opt-DMD');
set(h,{'LineWidth'},{0.5;0.5;2;2})
legend ('Hare','Lynx','DMD Hare','DMD Lynx')

saveas (figure(1),'opt_DMD.jpg')

%% BOP-DMD Application

Niter =     1000;              % #BOP iterations
Imax =      train;             % Max # of data
Nsubset =   floor(0.6*train);  % Subset of data
[Ndataset,~] = size(HL);       % x
e1 = e;                        % Initial condition
temp = [];

% BOP-DMD Iterations

for i=1:Niter
    
    ind = randperm(Imax, Nsubset);
    ind = sort(ind,2,"ascend");
    temp(:,:) = HL(:,ind);
    time_b = time(ind);
    [w1,e1,b1] = optdmd(temp,time_b,2,1,[],e1);  %Iterations of the BOP-DMD
    wx(:,:,i)=w1;
    ex(:,:,i)=e1;
    bx(:,:,i)=b1;

%     figure(2)
%     plot (real(ex(:,:,i)),imag(ex(:,:,i)),'o')
%     hold on
end

clear temp 

% Constraints
[~,c1] = find( abs(real(ex(1,:,:)))>0.03 );
[~,c2] = find( abs(real(ex(2,:,:)))>0.03 );
c = union(c1, c2);

ex(:,:,c)=[];
bx(:,:,c)=[];
wx(:,:,c)=[];

% Mean Value and variance of parameters

Wmean = [ mean(wx(1,1,:),3), mean(wx(1,2,:),3);mean(wx(2,1,:),3), mean(wx(2,2,:),3)];

clear i

time_prevision = [0:0.1:30];

for i=1:length(ex)
    l_prev(:,:,i) = exp(ex(:,:,i).*time_prevision);
    x_prevision(:,:,i) = Wmean(:,:)*l_prev(:,:,i).*bx(:,:,i);
end

X_prev_mean = [mean(x_prevision(1,:,:),3); mean(x_prevision(2,:,:),3)];
X_prev_varn = [var(x_prevision(1,:));var(x_prevision(2,:))];

H_prev(:,:) = X_prev_mean(1,:,:);
L_prev(:,:) = X_prev_mean(2,:,:);

timeprev = years(1)+time_prevision*2;

figure(3)
subplot(1,3,1)

plot (years,H,'ro--',years,L,'bo--','LineWidth',0.5)
hold on
plot(timeprev,H_prev,'r',timeprev,L_prev,'b','LineWidth',2)
fill([years(end)-6 years(end) years(end) years(end)-6],[0 0 max(H) max(H)],'c','FaceAlpha',0.1,'LineStyle','-')
grid on
xlim([1845 1903])
ylim([0 max(H)])
ylabel ('Population State');
xlabel ('Years');
title ('BOP-DMD');
legend ('Hare','Lynx','BOP-DMD Hare','BOP-DMD Lynx')




%% Time Delay Embedding

% Hankel Matrix
p = floor(0.7*train);   % New time interval
% m = length(HL)-p+1;     % Double data sets
m = train-p+1;
HL_tde = zeros(m,p);    % Data matrix
DT = 1;                 % Time step

for j=1:2:2*m
    for i=1:p
        temp(:,i) = HL_train(:,i+DT-1);
    end
    HL_tde(j:(j+1),:) = temp(:,:);
    DT = DT+1;    
end

[u_tde,s_tde,v_tde] = svd(HL_tde,'econ');

s_sum = sum(diag(s_tde));
diag_s = diag(s_tde);
n = length(diag_s);

figure(2)
plot(1:n, diag_s(1:end)/s_sum*100,'ro','LineWidth',2)
hold on
fill([0 5 5 0],[0 0 60 60],'w','FaceAlpha',0.001,'LineStyle','-')
grid on
ylabel('${\it} \frac{\sigma_n}{\sum_{i=1}^{N} \sigma_i} \% $','Interpreter','Latex')
xlabel('Singular Value')
xlim([1 20])

rank = 5;
time_tde = 1:p;

[w_tde,e_tde,b_tde,atilde_tde] = optdmd(HL_tde,time_tde,rank,2);
time_int=[0:0.1:30];

HLx_prev = w_tde*diag(b_tde)*exp(e_tde.*time_int);

timeint = years(1)+time_int*2;

figure(3)
subplot(1,3,2)

plot (years,H,'r--o',years,L,'b--o','LineWidth',0.5)
hold on
plot(timeint,abs(HLx_prev(1,:)),'r',timeint,abs(HLx_prev(2,:)),'b','LineWidth',2)
fill([years(end)-6 years(end) years(end) years(end)-6],[0 0 max(H) max(H)],'c','FaceAlpha',0.1,'LineStyle','-')
xlim([1845 1903])
ylim([0 max(H)])
grid on
ylabel ('Population State');
xlabel ('Years');
title ('TDE-DMD');
legend ('Hare','Lynx','TDE-DMD Hare','TDE-DMD Lynx')




figure(4)
subplot (2,1,1)
plot (years,H,'r--o',years,L,'b--o','LineWidth',0.5)
hold on
plot(timeint,abs(HLx_prev(1,:)),'r',timeint,abs(HLx_prev(2,:)),'b','LineWidth',2)
fill([years(end)-6 years(end) years(end) years(end)-6],[0 0 max(H) max(H)],'c','FaceAlpha',0.1,'LineStyle','-')
xlim([1845 1903])
ylim([0 max(H)])
grid on
ylabel ('Population State');
xlabel ('Years');
label = sprintf('TDE-DMD: rank = %d', rank);
title (label);
legend ('Hare','Lynx','TDE-DMD Hare','TDE-DMD Lynx')


rank = 7;
time_tde = 1:p;

[w_tde,e_tde,b_tde,atilde_tde] = optdmd(HL_tde,time_tde,rank,2);
time_int=[0:0.1:30];

HLx_prev = w_tde*diag(b_tde)*exp(e_tde.*time_int);

timeint = years(1)+time_int*2;

figure(4)
subplot (2,1,2)
plot (years,H,'r--o',years,L,'b--o','LineWidth',0.5)
hold on
plot(timeint,abs(HLx_prev(1,:)),'r',timeint,abs(HLx_prev(2,:)),'b','LineWidth',2)
fill([years(end)-6 years(end) years(end) years(end)-6],[0 0 max(H) max(H)],'c','FaceAlpha',0.1,'LineStyle','-')
xlim([1845 1903])
ylim([0 max(H)])
grid on
ylabel ('Population State');
xlabel ('Years');
label = sprintf('TDE-DMD: rank = %d', rank);
title (label);
legend ('Hare','Lynx','TDE-DMD Hare','TDE-DMD Lynx')

saveas (figure(4),'TDE_DMD_rank_variation.jpg')

%% Bagging- DMD Time Delay Embedding

Niter   = 500;
Imax    = p;
Nsubset = floor(0.9*p);
temp    = [];
w_n     = [];
e_n     = [];
b_n     = [];
e_bt    = e_tde;
temp    = [];

for i=1:Niter
    
    ind = randperm(Imax, Nsubset);
    ind = sort(ind,2,"ascend");
    temp(:,:) = HL_tde(:,ind);
    Time = time(ind);
    
    [w_bt,e_bt,b_bt] = optdmd(temp,Time,rank,1,[],e_tde); %Iterations of the BOP-DMD
    w_n(:,:,i) = w_bt;
    e_n(:,:,i) = e_bt;
    b_n(:,:,i) = b_bt;

%     figure(2)
%     plot (real(ex(:,:,i)),imag(ex(:,:,i)),'o')
%     hold on
end

clear temp 

[~,c1] = find( abs(real(e_n(1,:,:)))>0.03 );
[~,c2] = find( abs(real(e_n(2,:,:)))>0.03 );
c = union(c1, c2);
e_n(:,:,c)=[];
b_n(:,:,c)=[];
w_n(:,:,c)=[];

for i = 1:size(w_n,1)
    for j = 1:size(w_n,2)

        Wmean_bt(i,j,:) = mean(w_n(i,j,:),3);
    end
end

clear i

DeltaT = 1;

time_prevision_bt = [0:0.1:30];
x_prevision_bt = [];

for i=1:length(e_n)
    l_prev_bt(:,:,i) = exp(e_n(:,i).*time_prevision_bt);
    x_prevision_bt(:,:,i) = Wmean_bt*diag(b_n(:,:,i))*l_prev_bt(:,:,i);
end

X_prev_mean_bt = [mean(x_prevision_bt(1,:,:),3);mean(x_prevision_bt(2,:,:),3)];
X_prev_varn_bt = [var(x_prevision_bt(1,:));var(x_prevision_bt(2,:))];

H_prev_bt(:,:) = X_prev_mean_bt(1,:,:);
L_prev_bt(:,:) = X_prev_mean_bt(2,:,:);


timeintb = years(1)+time_prevision*2;

figure(3)
subplot(1,3,3)
plot (years,H,'r--o',years,L,'b--o','LineWidth',0.5)
hold on
plot(timeintb,H_prev_bt,'r',timeintb,L_prev_bt,'b','Linewidth',2)
fill([years(end)-6 years(end) years(end) years(end)-6],[0 0 max(H) max(H)],'c','FaceAlpha',0.1,'LineStyle','-')
xlim([1845 1903])
ylim([0 max(H)])
grid on
ylabel ('Population State');
xlabel ('Years');
title ('Bagging TDE-DMD');
legend ('Hare','Lynx','BTDE-DMD Hare','BTDE-DMD Lynx')

saveas (figure(3),'DMD_TDE_BOP_comparison.jpg')
