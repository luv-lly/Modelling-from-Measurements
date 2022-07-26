%% NN - KS Equations
clear all;close all;clc;

input=[];
output=[];

N = 4;            % States
N_iteration= 100; % #Iteration random conditions 


for j=1:N_iteration
    u0 = randn(N,1);
    [t,x,u] = ks_m(u0,N);

    input  = [input; u(1:end-1,:)];
    output = [output; u(2:end,:)];
    
end

figure(3)
pcolor(x,t,u), shading interp, colormap(hot)

%% Train NN
net = feedforwardnet([10 10 10]);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'logsig';
net.layers{3}.transferFcn = 'purelin';

net = train(net,input.',output.');

u1 = u(1,:).';
u_nn(1,:)=u1;
for i=2:length(t)

    utemp = net(u1);
    u_nn(i,:)=utemp.';
    u1=utemp;

end

%% Compare the results using different initial conditions
u_nc = randn(N,1);

[t_real,x_real,u_real] = ks_m(u_nc,N);

u_test1 = u_real(1,:).';
u_nn(1,:)=u_test1;

for jj=2:length(t_real)
    u_next = net(u_test1);
    u_nn(jj,:)=u_next.';
    u_test1=u_next;
end

 %%   
figure(2)
fig2 = figure(2);
subplot(1,2,1)
pcolor(x_real,t_real,u_real),shading interp, colormap(hot)
title('KS - Original data set')
xlabel('space (norm.)')
ylabel('time (norm.)')
subplot(1,2,2)
pcolor(x_real,t_real,u_nn),shading interp, colormap(hot)
title('KS - NN data')
xlabel('space (norm.)')
ylabel('time (norm.)')
figure(4)
pcolor(x_real,t_real,(u_real-u_nn)./u_real),shading interp, colormap(hot)
title('error')
colorbar


saveas(fig2,'KS_NN.jpg') 

%% Compare the results using different initial conditions
u_nc = randn(N,1);

[t_real,x_real,u_real] = ks_m(u_nc,N);

u_test1 = u_real(1,:).';
u_nn(1,:)=u_test1;

for jj=2:length(t_real)
    u_next = net(u_test1);
    u_nn(jj,:)=u_next.';
    u_test1=u_next;
end

 %%   
figure(3)
fig3 = figure(3);
subplot(1,2,1)
pcolor(x_real,t_real,u_real),shading interp, colormap(hot)
title({'KS - Original data set'; ' Different initial conditions'})
xlabel('space (norm.)')
ylabel('time (norm.)')
subplot(1,2,2)
pcolor(x_real,t_real,u_nn),shading interp, colormap(hot)
title({'KS - NN data'; ' Different initial conditions'})
xlabel('space (norm.)')
ylabel('time (norm.)')
figure(4)
pcolor(x_real,t_real,(u_real-u_nn)./u_real),shading interp, colormap(hot)
title('error')
colorbar


saveas(fig3,'KS_NN_NC.jpg') 