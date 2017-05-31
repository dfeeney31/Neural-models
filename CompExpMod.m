%% The most realistic set of data below
%As of 3/19, this creates a set of N neurons with a latent
%state . As a result, A is 2x2 identity with some noise and multiplied through
%the 500 iterations. C is Nx2 and encodes the effect of the latent state on
%each neuron (each row). C was estimated to create realistic discharge rates. 
%Xs is how X changes through the trial and includes
%Gaussian noise awgn. 
%z(t) = Cx(t). Discharge rate is exp(z(t)).
%Start by estimating the number of active motor units given a desired force
%Output of the model is stored in seq, which is run through PLDSExample

clear all
%% Make force trace 
GoalF = 5;
force(1:101) = 0:GoalF/100:GoalF;                      %Linear 1s ramp with 4 second hold of a target force
force(101:500) = force(101);

%% Parameters of model: X and A
Xs = zeros(500,2);
A = [1,0; 0, 1];                        %This is the linear approximation of the dynamics and should always be dxd where d is dimensionality of latent state
Xi = [1 ; 1];                           %This is the initial latent state
Xs(1,:) = A * Xi;                       %How the latent state evolves through the trial

%% Estimate the numbers of active motor units for the 5 second trial. Based on Fuglevand, 1993 on the excitation level a unit begins
%Discharging action potentials. The force in this case is a proxy for the
%synaptic current

RTE = zeros(1,120); %Preallocating the vector of recruitment excitation levels. Sorry for the for loop
a = log(60)/120;    %Parameters of exponential relation. Last unit recruited at 60% MVC
for ff = 1:120
    RTE(ff) = exp(a*ff);
end

%Need to calcualte, at each instant, how far above recruitment threshold a
%unit is
AboveRTE = zeros(120,500);                  %Preallocate matrix where each row is a neuron and length of 500 ms.
for tt = 1:500
    AboveRTE(:,tt) = force(tt) - RTE;          %Find the difference between recruitment threshold and current force output
end

ActiveUnits =zeros(1,120);          %Preallocate how many units will eventually be 'turned on' by this force

for yy = 1:500
    nums = find(AboveRTE(:,yy)>0);                %Find all entries above zero
    ActiveUnits(yy) = length(nums);     %Find the index where RTE is greater than the force.This is approx # motor units active 
end
AboveRTE(AboveRTE<0)=0; %If units are not above RTE during the trial, they should be disregarded
neurons = ActiveUnits(500); %Total number of units discharging action potentials during the trial

%% Set the yDim to the number of motor activated (from the previous step)
Trials = 1; yDim = ActiveUnits(500); T = 500;
seq = struct();

%% The following will need to adjust the firing rates for the N active motor units using x, A, and u
for i = 2:500
    Xs(i,:) = A * Xs(i-1,:)';  %Evolution of latent state. Currently a boring identity matrix
end
Xs = awgn(Xs,35);                %Adding Gaussian noise to evolution, which models synaptic fluctations 
 
C = zeros(yDim,2);       %making driving matrix C be n rows (one for each neuron) and two columns (one for each latent input)
for mm = 1:yDim
    C(mm,1) = 2 + 0.00541*exp((0.049)*mm);        %Set the linear equation for how C changes for each neuron for input conductance
    C(mm,2) = 0.5 + 0.001666 * exp((0.048)*mm);     %Set linear equaiton for C for AHP duration
end


LinMult = AboveRTE ./5;                         %Create column array of force differences to model sigmoid relation between input force and discharge rate
LinMult = tanh(LinMult) * 0.5 + 0.5;             
LinMult = LinMult(1:ActiveUnits(500),:);            %Get to correct size

FiringRate = zeros(ActiveUnits(500),500);           %preallocate
AboveRTE = AboveRTE(1:ActiveUnits(500),1:500);
NoiseMat = ones(ActiveUnits(500),1); 
NoiseMat = NoiseMat .* (AboveRTE(:,500).* 5);
NoiseMat(NoiseMat<0) = 0;


%% Use the noise plut the discharge rates and send through exponential.
for kk = 1:500
    FiringRate(:,kk) = (C * Xs(kk,:)').* LinMult(:,kk); %Average firing rate for each neuron in column array 
end

for ooo = 1:ActiveUnits(500)
    FiringRate(ooo,:) =  awgn(FiringRate(ooo,:), NoiseMat(ooo)); %Add noise to the firing rate
end

FiringRate = exp(FiringRate); %Send through exponential

MultFR = FiringRate;

for rr = 1:ActiveUnits(500)
    for ss = 1:length(AboveRTE)
        if AboveRTE(rr,ss) == 0;
            MultFR(rr,ss) = 0;
        end
    end
end
AvgFiringRate = mean(MultFR(:,300:500),2);
%% Find the lambda value for PoissonRnd.
rng(123);
DummyVariable = zeros(ActiveUnits(500),5); %This is going to be the 5 seconds drawn from the Poisson Distribution
for JA = 1:ActiveUnits(500)
    DummyVariable(JA,:) = poissrnd(AvgFiringRate(JA),1,5);
end

Try = zeros(67,500);
for dan = 1:ActiveUnits(500)
   Try(dan,1:100) = binornd(ones(1,100),DummyVariable(dan,1)/100); 
   Try(dan,101:200) = binornd(ones(1,100),DummyVariable(dan,2)/100);
   Try(dan,201:300) = binornd(ones(1,100),DummyVariable(dan,3)/100);
   Try(dan,301:400) = binornd(ones(1,100),DummyVariable(dan,4)/100);
   Try(dan,401:500) = binornd(ones(1,100),DummyVariable(dan,2)/100);
end


seq.y = Try;
seq.T = T;
seq.yr = MultFR;
seq.x = Xs;
%seq.u = ones(2,500).*0.5; Unsure the best inputs for u

NumPlots = floor(ActiveUnits(500)/10); 
figure(1)
for aa = 1:NumPlots
    subplot(3,5,aa);
    plot(1:500, MultFR(aa*10,:),'*');
    xlabel('bin');
    ylabel('Firing Rate Hz');
    axis([0,500,0,20]);
end

%% print DR and CV values for time 0.5 - 1 s and 3.5 to 4 s
Multiplier = floor(ActiveUnits(500)/10);

ISI = 1./MultFR;  %Calculate the interspike interval at each time instant
ISI2 = ISI;
ISI(ISI==inf)=0;
ISI2(ISI2==inf)=0;
ISI2 = 10^3 .*ISI2;               %Convert to ms from seconds

%%% Time vector voodoo %%%
rounded_ISI = floor(ISI2);
rounded_ISI(isnan(rounded_ISI)) = 0; 
rounded2 = rounded_ISI(:,25:end);
for dd = 1:length(AboveRTE(:,500))
    discharge_index(dd,:) = cumsum(rounded_ISI(dd,:));
end


SDlow = zeros(1,Multiplier);
AVlow = zeros(1,Multiplier);
CVlow = zeros(1,Multiplier);
SDhi = zeros(1,Multiplier);
AVhi = zeros(1,Multiplier);
CVhi = zeros(1,Multiplier);

SDlow(1) = (std(ISI(1,50:100)))^2;          %Calculating SD^2 for ISI
AVlow(1) = (mean(ISI(1,50:100)))^3;         %Calculating avg^3 for ISI
CVlow(1) = sqrt((SDlow(1)/AVlow(1)))*100;
SDhi(1) = (std(ISI(1,350:400)))^2;
AVhi(1) = (mean(ISI(1,350:400)))^3;
CVhi(1) = sqrt((SDhi(1)/AVhi(1)))*100;
DRlow(1) = mean(MultFR(1,50:100));
DRhi(1) = mean(MultFR(1,350:400));
outcomes(1,:) = [AVlow(1), CVlow(1), AVhi(1),CVhi(1)];

for rr = 1:Multiplier
    SDlow(rr) = (std(ISI(rr*10,101:201)))^2;
    AVlow(rr) = (mean(ISI(rr*10,101:201)))^3;
    CVlow(rr) = sqrt((SDlow(rr)/AVlow(rr)))*100;
    SDhi(rr) = (std(ISI(rr*10,350:400)))^2;
    AVhi(rr) = (mean(ISI(rr*10,350:400)))^3;
    CVhi(rr) = sqrt((SDhi(rr)/AVhi(rr)))*100;
    DRlow(rr) = mean(MultFR(rr*10,101:201));
    DRhi(rr) = mean(MultFR(rr*10,350:400));
    outcomes(rr,:) = [DRlow(rr), CVlow(rr), DRhi(rr),CVhi(rr)];
end

Stats = {'DRi' 'CVi' 'DRf' 'CVf'};                  %Stats outputs a cell array with Discharge rate and CV initially (0.5-1 s) and at end (3.5- 4 s)
Stats(2:Multiplier+1,1) = num2cell(DRlow);          %(It's a cell array so each column can be labeled)
Stats(2:Multiplier+1,2) = num2cell(CVlow);
Stats(2:Multiplier+1,3) = num2cell(DRhi);
Stats(2:Multiplier+1,4) = num2cell(CVhi);

%% Isometric force production model from motor unit firing times

P = zeros(1,ActiveUnits(500));   %Preallocate vector of peak twitch force for each motor unit
b = log(100)/120;            %Constant multiplier for setting up force multiplier
for nn = 1:ActiveUnits(500)
    P(nn) = exp(b*nn);
end

TwitchTimes = zeros(1,ActiveUnits(500));  %Preallocate twitch contraction times
c = logb(100,3);
d = 1/c;
for ccc = 1:ActiveUnits(500)
    TwitchTimes(ccc) = 90*((1/P(ccc))^d);
end

NormTi = zeros(ActiveUnits(500),length(Xs));        
for hhh = 1:length(Xs)
    NormTi(:,hhh) = (TwitchTimes')./ISI2(:,hhh);
end
NormTi(NormTi==inf)=NaN;

force_mat = zeros(ActiveUnits(500),length(Xs));  %Preallocate the force for each motor unit at each time point
PoverT = P./TwitchTimes;
for fff = 1:length(Xs)
    force_mat(:,fff) = PoverT.* exp(1-(1/TwitchTimes')); %Gives the force output with no variability
end

S = -2.*((NormTi).^3);        %Makes the sigmoid to determine gain
S = 1 - exp(NormTi);
G = abs(S./NormTi);
G(G==inf)=0;

for rrrr = 1:ActiveUnits(500)
    for cccc = 1:length(Xs)
        if NormTi(rrrr,cccc) > 0.4
            force_mat(rrrr,cccc) = force_mat(rrrr,cccc) * G(rrrr,cccc);
        end
    end
end

for hello = 1:ActiveUnits(500)
    for bye = 1:length(Xs)
        if AboveRTE(hello,bye) == 0;
           force_mat(hello,bye) = 0;
        end
    end
end

Total_Force = sum(force_mat);
figure(2);
plot(Total_Force);          %Scaling up to arbitrary units
% ylim([0,10]);
ylabel('Force (au)');
xlabel('Time (ms)');
mean_force = mean(Total_Force(301:500));
SD_force = std(Total_Force(301:500));
CV_force = (SD_force./mean_force) * 100

%% Find a spectral density estiamte. 

Fout = Total_Force;
FS = 100;
figure(3);
pmtm(Fout,4,length(Fout),FS);
OutputDR = mean(MultFR(:,300:500),2);

%% Interval histograms

ISI2(ISI2 == 0) = NaN; 
ISI2(ISI2 > 200) = NaN;
figure(4)
histogram(ISI2(47,:))

Kurtosis = kurtosis(ISI2(47,300:500));
Skewness = skewness(ISI2(47,300:500));

Kurt_Value = nanmean(Kurtosis);
Skew_Value = nanmean(Skewness);
%dlmwrite('Model_Output2.csv',AvgFiringRate','-append');

figure(5);
L = logical(Try);
plotSpikeRaster(L,'PlotType','Scatter');

%% Raster 
rng(123);

lambda = zeros(neurons, length(force));
for j = 1:length(force)
    for mmmm = 1:neurons
        lambda(mmmm,j) = MultFR(mmmm,j);    
    end 
end
   
time = floor(exprnd(lambda)); 
Spike_vec = zeros(neurons,5000);
for lm = 1:neurons
    for jj = 10:length(force)
        if time(lm,jj) > 0
            Spike_vec(lm,(jj*10)+time(lm,jj)) = 1;
        else
            Spike_vec(lm,(jj*10)+time(lm,jj)) = 0;
        end
    end
end

newL = floor(length(Spike_vec)/10)
for nn = 1:neurons
    for d = 0:(newL-1)
        dsSub(nn,d+1) = sum(Spike_vec(nn,(d*10)+1:(d*10)+10));
    end
end

seq.y = dsSub
seq.T = length(dsSub);
seq.yr = MultFR;
%seq.u = ones(2,500).*0.5; Unsure the best inputs for u

figure;
plotSpikeRaster(logical(dsSub),'PlotType','Scatter');
set(gcf,'color','white')