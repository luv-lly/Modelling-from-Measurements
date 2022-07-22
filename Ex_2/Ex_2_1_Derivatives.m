
%% 2nd Exercise
% Neural Network KS equations
clear all; close all; clc;

% Space/time domain
load('kuramoto_sivashinsky.mat')
x = xsave;
t = tsave;
u = usave;

% Normalization of dataset between 0 and 1

max_v = max(max(usave));
min_v = min(min(usave));

scale_down = @(x) (x - min_v) / (max_v - min_v);
scale_up = @(x) x * (max_v - min_v) +  min_v;

usave = scale_down(usave);

% Original DATA plot
%  figure (3)
%  pcolor(xsave,tsave,usave),shading interp, colormap(hot)

ud  = [];
ud2 = [];
ud3 = [];
ud4 = [];

xlast = ones([1024,1]);
usave1 = (cat (2,usave,usave));
usave1_t = usave1.';
xsave1 = (cat (1,xsave, xlast+xsave));

for i = 1: length (t)

    temp1 = diff(usave1_t(:,i))./diff(xsave1(:,1));
    ud = [ud, temp1];
    temp2 = diff(ud(:,i))./diff(xsave1(2:end,1));
    ud2 =[ud2, temp2];
    temp3 = diff(ud2(:,i))./diff(xsave1(3:end,1));
    ud3 =[ud3, temp3];
    temp4 = diff(ud3(:,i))./diff(xsave1(4:end,1));
    ud4 =[ud4, temp4];
end

ud = ud.';
ud2 = ud2.';
ud3 = ud3.';
ud4 = ud4.';

input  = [];
output = [];

for i = 1:length(xsave')
    
    ux =usave1(1:70,i+4);
    uy =usave1(2:71,i);
    input  = cat (1, input , [usave1(1:70,i+4) ud(1:70,i+4) ud2(1:70,i+4) ud3(1:70,i+4) ud4(1:70,i+4)]);
    output = cat (1, output, usave1(2:71,i));

end


% Train the NN
HiddenLayer = [15 10 5];
TrainFNC = 'trainlm';

net = fitnet(HiddenLayer, TrainFNC);

 net.divideParam.trainRatio = 85/100;
 net.divideParam.valRatio = 20/100;
 net.divideParam.testRatio = 20/100;

 net.output.processFcns = {'removeconstantrows','mapminmax'};
 net.divideMode = 'sample';

[net, tr] = train(net,input.', output.');

% Test the NN

y = net ( input.' );
y = reshape(y',[70,1024]);

y = scale_up(y);


figure (1)
subplot(1,2,1)
pcolor(xsave,tsave(2:end),usave(1:70,:)), shading interp, colormap(hot)
ylabel ('time');
xlabel ('x');
title ('KS data set');

subplot(1,2,2)
pcolor(xsave,tsave(2:end),y), shading interp, colormap(hot)
ylabel ('time');
xlabel ('x');
title ('KS NN');

saveas(figure(1),'KS_NN_derivatives.jpg') 

figure(2)
plot (xsave, y(9,:))
hold on
plot (xsave, u(10,:))

ylabel ('x');
xlabel ('u');
title ('KS evolution');

saveas(figure(1),'KS_evolution_derivatives.jpg') 



