%% NN - Lorenz equations
clear all; close all; clc;

% Standard coefficients
b   = 8/3;
sig = 10;

% Time series
T = 10;
dt = 0.01;
t = 0:dt:T;

%Lorenz equations in the different conditions

rho_train =[10 28 35];
rho_test = [17 40];

Lrnz_1 = @(t,x) ([sig*(x(2)-x(1)); rho_train(1)*x(1)-x(1)*x(3)-x(2); x(1)*x(2)-b*x(3)]);
Lrnz_2 = @(t,x) ([sig*(x(2)-x(1)); rho_train(2)*x(1)-x(1)*x(3)-x(2); x(1)*x(2)-b*x(3)]);
Lrnz_3 = @(t,x) ([sig*(x(2)-x(1)); rho_train(3)*x(1)-x(1)*x(3)-x(2); x(1)*x(2)-b*x(3)]);

Lrnz_4 = @(t,x) ([sig*(x(2)-x(1)); rho_test(1)*x(1)-x(1)*x(3)-x(2); x(1)*x(2)-b*x(3)]);
Lrnz_5 = @(t,x) ([sig*(x(2)-x(1)); rho_test(2)*x(1)-x(1)*x(3)-x(2); x(1)*x(2)-b*x(3)]);


ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
input  = [];
output = [];

% NN Training phase

% Data generation
for j=1:50  

x0=30*(rand(3,1)-0.5);
    [t,y] = ode45(Lrnz_1,t,x0);
    rho = ((y(2:end,2)-y(1:end-1,2))/dt+y(1:end-1,1).*y(1:end-1,3)+y(1:end-1,2))./y(1:end-1,1);
    input=[input; y(1:end-1,:) rho.*y(1:end-1,1)];
    output=[output; y(2:end,:)];
    plot3(y(:,1),y(:,2),y(:,3)), hold on
    plot3(x0(1),x0(2),x0(3),'ro')
    
    [t,y] = ode45(Lrnz_2,t,x0);
    rho = ((y(2:end,2)-y(1:end-1,2))/dt+y(1:end-1,1).*y(1:end-1,3)+y(1:end-1,2))./y(1:end-1,1);
    input=[input; y(1:end-1,:) rho.*y(1:end-1,1)];
    output=[output; y(2:end,:)];
    plot3(y(:,1),y(:,2),y(:,3)), hold on
    plot3(x0(1),x0(2),x0(3),'ro')
    
    [t,y] = ode45(Lrnz_3,t,x0);
    rho = ((y(2:end,2)-y(1:end-1,2))/dt+y(1:end-1,1).*y(1:end-1,3)+y(1:end-1,2))./y(1:end-1,1);
    input=[input; y(1:end-1,:) rho.*y(1:end-1,1)];
    output=[output; y(2:end,:)];
    plot3(y(:,1),y(:,2),y(:,3)), hold on
    plot3(x0(1),x0(2),x0(3),'ro')
         
end
   
%% NN Train
trainFcn = 'trainlm';
net = fitnet([10 10],trainFcn);

net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample

net.divideParam.trainRatio = 80/100;
net.divideParam.valRatio = 25/100;
net.divideParam.testRatio = 15/100;

% Train the NN
[net,tr] = train(net, input.',output.');

%% NN Test phase

%% Test for rho = 17
x0=30*(rand(3,1)-0.5);
[t,y] = ode45(Lrnz_4,t,x0);

rho_nn(1)=((y(2,2)-y(1,2))/dt+y(1,1).*y(1,3)+y(1,2))./y(1,1); %rho initial condition
for i = 2:length(t)
    % Recursive NN
    y0= net([x0(1:3); rho_nn(i-1)*x0(1)]);
    y_nn(i,1:3) = y0(1:3).';
    rho_nn(i) = ((y(i,2)-y(i-1,2))/dt+y(i-1,1).*y(i-1,3)+y(i-1,2))./y(i-1,1);
    x0=y0;

end

figure(1)

subplot(3,2,1), plot(t,y(:,1),'b',t,y_nn(:,1),'r','Linewidth',2)
title ('Lorenz equation - \rho = 17');
xlabel ('time');
ylabel ('x');
legend ('Data','Forecast');
grid on
subplot(3,2,3), plot(t,y(:,2),'b',t,y_nn(:,2),'r','Linewidth',2)
xlabel ('time');
ylabel ('y');
legend ('Data','Forecast');
grid on
subplot(3,2,5), plot(t,y(:,3),'b',t,y_nn(:,3),'r','Linewidth',2)
xlabel ('time');
ylabel ('z');
legend ('Data','Forecast');
grid on



%% Test for rho = 40

x0=30*(rand(3,1)-0.5);
[t,y] = ode45(Lrnz_5,t,x0);

rho_nn(1)=((y(2,2)-y(1,2))/dt+y(1,1).*y(1,3)+y(1,2))./y(1,1); %rho initial condition
for i=2:length(t)
    % Recursive NN
    y0=net([x0(1:3); rho_nn(i-1)*x0(1)]);
    y_nn(i,1:3)=y0(1:3).';
    rho_nn(i)=((y(i,2)-y(i-1,2))/dt+y(i-1,1).*y(i-1,3)+y(i-1,2))./y(i-1,1);
    x0=y0;
end

figure(1)
subplot(3,2,2), plot(t,y(:,1),'b',t,y_nn(:,1),'r','Linewidth',2)
title ('Lorenz equation - \rho = 40');
xlabel ('time');
ylabel ('x');
legend ('Data','Forecast');
grid on
subplot(3,2,4), plot(t,y(:,2),'b',t,y_nn(:,2),'r','Linewidth',2)
xlabel ('time');
ylabel ('y');
legend ('Data','Forecast');
grid on
subplot(3,2,6), plot(t,y(:,3),'b',t,y_nn(:,3),'r','Linewidth',2)
xlabel ('time');
ylabel ('z');
legend ('Data','Forecast');
grid on


saveas(figure(1),'Lorenz_forecast.jpg')

