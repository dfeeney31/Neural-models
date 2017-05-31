%% Create 2 LIF neurons with a sinusoidal input current where neuron 1 inhibits neuron 2 with each discharge of an AP
% Created Dan Feeney Jan, 2017 based on Gao 2007 and others. Units are AU
% unless otherwise noted. 

clear
% input current
f = 3; %frequency of stimulation

% capacitance and leak resistance
C = 1; % nF
R = 40; % mohms

% I & F implementation dV/dt = - V/RC + I/C
% Using h = 1 ms step size, Euler method

y = [0:1/1000:1];
A = 1;
f = 3;
I = A*sin(2*pi*f*y) + 0.5;
plot(y,I)

rng(10);
V = 0;
V2 = 0;
tstop = 1000;
I = awgn(I,45);

abs_ref = 50; % AHP duration 
ref = 0; % absolute refractory period counter
ref2 = 0; %counter neuron 2
V_trace = []; % voltage trace for plotting
V_trace2 = []; %Second neuron
V_th = 10; % spike threshold
V_peak = 50; %Spike height

%% Run through time for neuron 1, then for neuron 2
for t = 1:tstop
  
   if ~ref
     V = V - (V/(R*C)) + (I(t)/C);
     V2 = V2 - (V2/(R*C)) + (I(t)/C); 
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

ind = (V_trace == 50); %Return the indices where neuron 1 discharged an AP.
%% Run for neuron 2 that will receive (-)ive current for each discharge of neuron 1
%Set the negative current to inject into neuron 2 at each index of neuron 1 discharging

I(end) = [];
I2 = I + ind;
for i = 1:length(I)
   if (ind(i) == 1)
       I2(i:i+10) = -0.5:0.05:0;
   end
end
I2 = awgn(I2,10);
I2 = I2 + 0.5;

for t2 = 1:tstop
   if ~ref2
     V2 = V2 - (V2/(R*C)) + (I2(t2)/C); 
   else
     ref2 = ref2 - 1;
     V2 = 0.2*V_th; % reset voltage
   end
   
   if (V2 > V_th)
     V2 = V_peak;  % AP
     ref2 = abs_ref; % set refractory counter
   end

   V_trace2 = [V_trace2 V2];
end

  figure(1)
  plot(V_trace)
  figure(2)
  plot(V_trace2)
  
  %% print out some summary statistics
  num_spikes = sum(V_trace==V_peak);
  spikes2 = sum(V_trace2==V_peak);
  DR = num_spikes/(tstop/1000);
  DR2 = spikes2/(tstop/1000);
  pks = length(findpeaks(I));
  
  %% Make raster of discharges
  V1 = (V_trace==50);
  V2 = (V_trace2==50);
  output = cat(1, V1, V2);
  figure(3)
  plotSpikeRaster(output, 'PlotType', 'scatter')
  
seq = struct();
seq.y = double(output);
seq.T = 1000;
seq.u = ones(2,1000).*0.5;

%% Below is after the latent state trajectories are calcualted from PLDSExample
mscohere(params.model.Xpca(1,:),downsample(I,10),[],[],[],10)
xlim([0 5])
%mscohere(params.model.Xpca(2,:),downsample(clI,10),[],[],[],10)
%plotSpikeRaster(output,'plottype','scatter')
%plot(params.model.Xpca(1,:),'LineWidth',2,'Color',[0 0 0])
% plot(params.model.Xpca(2,:),'LineWidth',2,'Color',[0 0 0])
grid off
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);

plotyy(1:100,params.model.Xpca(1,:),1:100,downsample(I,10))
plot(1:100,params.model.Xpca(1,:))
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);

%% Separately compre frequency content of signals
Fs = 10
[P1,f1] = periodogram(params.model.Xpca(1,:),[],[],Fs,'power');
[P2,f2] = periodogram(I,[],[],1000,'power');

subplot(2,1,1)
plot(f1,P1,'k')
grid
ylabel('P_1')
title('Power Spectrum')
grid off
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);

subplot(2,1,2)
plot(f2,P2,'r')
grid
ylabel('P_2')
xlabel('Frequency (Hz)')
xlim([0 5])
grid off
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);


%% Low pass filter CST
CST = double(sum(output))
fc = 10;  %cutoff freq
fn = 1000/2; %For kinetic data at 100Hz
[b,a]=butter(6,fc/fn,'low');
cst_filt = filtfilt(b,a,CST);
plot(cst_filt)

[s c ph ci phi]= cmtm(cst_filt,I(1:1000),0.01,8,0,0,1)
plot(s,c,'k','linewidth',2);
xlim([1 10])
xlabel('Frequency')
ylabel('Coherence')
box off
set(gcf,'color','w')
set(gca,'TickDir', 'out','LineWidth',2,'TickDir','out','FontSize',16);