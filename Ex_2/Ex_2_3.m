clear all; close all; clc;
%% Reaction - diffusion equations

load ('reaction_diffusion_big.mat')

% Two-component reaction–diffusion equations (u,v)
% t = time sequence
% x = state sequence
% y = state sequence
% u = time-state evaluation
% v = time-state evaluation


%% Dimension reduction
nu = size(u,2);         % Number of columns of u
f_ux = zeros(0);        % Initialization of the flat matrix

for i = 1:length(t)

    f_ux_i = reshape(u(:,:,i),[1,nu^2]);
    f_ux = [f_ux;f_ux_i];

end
f_ux = f_ux';
[u_rd,s_rd,v_rd] = svd(f_ux,'econ');

% Plot Sigma
%C Percentage of variance of each mode

s_sum = sum(diag(s_rd));
diag_s = diag(s_rd);
n = length(diag_s);

figure(1)
plot(1:n, diag_s(1:end)/s_sum*100,'ro','LineWidth',2)
hold on
fill([0 10 10 0],[0 0 60 60],'w','FaceAlpha',0.001,'LineStyle','-')
grid on
ylabel('${\it} \frac{\sigma_n}{\sum_{i=1}^{N} \sigma_i} \% $','Interpreter','Latex')
xlabel('Singular Value')
xlim([1 20])

clear s m n

% Select rank
r = 10;

% Truncation up to the rank
u_tr = u_rd(:,1:r); 
s_tr = s_rd(1:r,1:r); 
v_tr = v_rd(:,1:r);

t_ind = round(0.95*length(t));
t_n =t(1:t_ind);
reduced_states = zeros(length(t_n), r);

for i = 1:length(t_n)
    f_state = f_ux(:,i);
    reduced_state_i = u_tr\f_state;
    reduced_states(i,:) = reduced_state_i;

end


%% Train NN

input   = reduced_states(1:end-1,:)';
output  = reduced_states(2:end,:)';

hiddenLayerSize = [10, 7, 5];
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.
net = fitnet(hiddenLayerSize,trainFcn);

 net.layers{1}.transferFcn = 'logsig';
 net.layers{2}.transferFcn = 'logsig';
 net.layers{3}.transferFcn = 'purelin';

% Choose Input and Output Pre/Post-Processing Functions
% For a list of all processing functions type: help nnprocess
net.input.processFcns = {'removeconstantrows','mapminmax'};
net.output.processFcns = {'removeconstantrows','mapminmax'};

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Train the Neural Network
[net,tr] = train(net,input,output);

%% Test NN
y_nn = net(input);

f_nn = [];

for i = 1:length(t_n)-1
    forecast_i = u_tr*y_nn(:,i);
    forecast_im = reshape(forecast_i, [nu,nu]);
    f_nn(:,:,i) = forecast_im;
    
end 

% compare frames
nframe = 123;

figure(2)
subplot(1,2,1)
pcolor(x,y,u(:,:,nframe)); shading interp; colormap(hot)
title(['Data set - time frame = ', num2str(nframe)])
xlabel('x')
ylabel('y')
subplot(1,2,2)
pcolor(x,y,f_nn(:,:,nframe)); shading interp; colormap(hot)
title(['NN forecast - time frame =', num2str(nframe)])
xlabel('x')
ylabel('y')

saveas(figure(2),'RD_NN.jpg')

%% Different Frame

T = t(11:end);
reduced_states = zeros(length(T), r);

for i = 1:length(T)
    f_state = f_ux(:,i);
    reduced_state_i = u_tr\f_state;
    reduced_states(i,:) = reduced_state_i;

end
input   = reduced_states(1:end-1,:)';
y_nn = net(input);

f_nn = [];

for i = 1:length(T)-1
    forecast_i = u_tr*y_nn(:,i);
    forecast_im = reshape(forecast_i, [nu,nu]);
    f_nn(:,:,i) = forecast_im;
    
end 

nframe = 192;

figure(3)
subplot(1,2,1)
pcolor(x,y,u(:,:,nframe)); shading interp; colormap(hot)
title(['Data set - time frame = ', num2str(nframe)])
xlabel('x')
ylabel('y')
subplot(1,2,2)
pcolor(x,y,f_nn(:,:,182)); shading interp; colormap(hot)
title(['NN forecast - time frame =', num2str(nframe)])
xlabel('x')
ylabel('y')

saveas(figure(3),'RD_NN_DifferentFrame.jpg')
