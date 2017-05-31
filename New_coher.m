%% Create 10 Integrate-Fire neurons with a sinusoidal input current with slightly different characteristics
% Created Dan Feeney Jan, 2017

clear
rng(10)
% input current
f = 4; %Frequency of waveform
f2 = 3;
% capacitance and leak resistance
C = [1:0.1:1.9]; % nF
R = [35:1:44]; % mohms

% I & F implementation dV/dt = - V/RC + I/C
% Using h = 1 ms step size, Euler method

y = [0:1/1000:1];
A = 1;
I = A*sin(f*2*pi*y) + 0.5; I = awgn(I,25);
I2 = A*sin(f2*2*pi*y) + 0.5; I2 = awgn(I2,25);
plot(y,I,'color',[0 0 0],'linewidth',2); hold on; plot(y,I2,'--k','linewidth',2)
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);
box off
set(gcf,'color','w')

V = 0;
V2 = 0;
tstop = 1000;

abs_ref = 80; % absolute refractory period 
ref = 0; % absolute refractory period counter
ref2 = 0; %counter neuron 2
V_trace = []; % voltage trace for plotting
V_th = 10; % spike threshold
V_peak = 50; %Spike height

%% Run through time for neuron 1, then for neuron 2
for neurons = 1:10
    
    for t = 1:tstop
        
        if ~ref
            V = V - (V/(R(neurons)*C(neurons))) + (I(t)/C(neurons));
            V2 = V2 - (V2/(R(neurons)*C(neurons))) + (I(t)/C(neurons));
        else
            ref = ref - 1;
            V = 0.2*V_th; % reset voltage
        end
        
        if (V > V_th)
            V = V_peak;  % Discharge AP
            ref = abs_ref; % set refractory counter
        end
        
        V_trace = [V_trace V];
        
    end
end

V_out = zeros(10,1000);
for B = 1:10
    j = (B-1)*1000;
    V_out(B,:) = V_trace(j+1:1000*B); 
end
%Plot the discharge times of the 10 units
for jj = 1:10
plot(V_out(jj,:))
hold on
end
hold off

ind = zeros(10,1000); %Preallocate binary vector of discharge instances
for l = 1:10
    ind(l,:) = (V_out(l,:) == 50); %Return the indices where neuron 1 discharged an AP.    
end
% MarkerFormat = struct()
%               MarkerFormat.Color = [0 0 0];
% plotSpikeRaster(logical(ind), 'plottype','scatter','markerformat',MarkerFormat)
% box off
% set(gcf,'color','w')
% set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);

%% Make a cumulative spike train for groups of 5 units
CST1 = sum(ind(1:5,:),1);
CST2 = sum(ind(6:10,:),1);
CSTtot = sum(ind,1);

plot(CST1); hold on; plot(CST2); plot(I);
%% Low pass filter the CST
fc = 10;  %cutoff freq
fn = 1000/2; %For kinetic data at 100Hz
[b,a]=butter(6,fc/fn,'low');
new_cst = filtfilt(b,a,CSTtot);
plot(new_cst)

[s4 c4 ph ci phi]= cmtm(new_cst,I(1:1000),0.01,8,0,0,1)
plot(s4,c4); xlim([0 10])
save('coher.mat','s4','c4','-append')
%plot(xs); hold on; plot(b);

[s c ph ci phi] = cmtm(xs,b,0.01,8,0,0,1);
s2 = s; c2= c;
save('coher.mat','s2','c2','-append')
xs2 = params.model.Xpca(1,:);
save('xs.mat','xs2','-append')
load('xs.mat')
  %% Make into sequence for PLDS model & plot
seq = struct();
seq.y = double(ind);
seq.T = 1000;
%seq.u = ones(2,1000).*0.5;
% 
plot(params.model.Xpca(1,:),'LineWidth',2,'Color',[0 0 0]); hold on
plot(xs,'k','LineWidth',2); hold on; plot(xs2,'--')
%hold on
%plot(downsample(I,10),'LineWidth',2,'Color',[0.25 0.25 0.25])
%plot(CST1,'LineWidth',2,'Color',[0 0 0])
% hold on
%plot(CST2,'LineWidth',2,'Color',[0 0 0])
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);
xlim([0 100])

mscohere(params.model.Xpca,downsample(I(1:1000),10),[],[],[],100)
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);
grid off
box off

ci(1:10) = 0.41;
plot(s3,c3,'k','linewidth',2)
hold on
plot(s4,c4,'--k','linewidth',2)
plot(ci,'--k','linewidth',2)
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);
xlim([1 10])

plot(I, 'k')
hold on
plot(I2,'k-')
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);
xlim([0 1000])
