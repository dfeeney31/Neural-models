%% The most realistic set of data below
%As of 12/16, this creates a set of N neurons with a latent
%state that is 2 dimensional. As a result, A is 2x2 identity with some noise and multiplied through
%the 500 iterations. C is Nx2 and encodes the effect of the latent state on
%each neuron (each row). C was estimated to create realistic discharge rates. 
%Xs is how X changes through the trial and includes
%Gaussian noise awgn. 
%z(t) = Cx(t). Discharge rate is exp(z(t)).
%Start by estimating the number of active motor units given a desired force

clear all
%% Make force trace 
GoalF = 7;
force(1:101) = 0:GoalF/100:GoalF;                      %Linear 1s ramp with 4 second hold of a target force
force(101:500) = force(101);

%% Parameters of model: X and A
Xs = zeros(500,2);
A = [1,0; 0, 1];                        %This is the linear approximation of the dynamics and should always be dxd where d is dimensionality of latent state
Xi = [1 ; 1];                           %This is the initial latent state
Xs(1,:) = A * Xi;                       %How the latent state evolves through the trial

%% Estimate the numbers of active motor units for the 5 second trial. Based on Fuglevand, 1993 on the excitation level a unit begins
%Discharging action potentials. The force in this case is a proxy for the
%synaptic current. Once it passes a percentage of max contraction (MVC),
%the unit is 'turned on'

RTE = zeros(1,120); %Preallocating the vector of recruitment excitation levels. Sorry for the for loop
a = log(60)/120;    %Parameters of exponential relation. Last unit recruited at 60% MVC in this case.
for ff = 1:120
    RTE(ff) = exp(a*ff);
end

%Need to calcualte, at each instant, how far above recruitment threshold a
%each motor unit is. 
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

%% Set the yDim to the number of motor activated (from the previous step). Seq is fed into decomposition algorithm.
Trials = 1; yDim = ActiveUnits(500); T = 500;
seq = struct();

%% The following will need to adjust the firing rates for the N active motor units using x, A, and u
for i = 2:500
    Xs(i,:) = A * Xs(i-1,:)';                   %Evolution of latent state. Currently a boring identity matrix
end
Xs = awgn(Xs,35);                               %Adding Gaussian noise to evolution, which models synaptic fluctations 
 
C = zeros(yDim,2);                              %making driving matrix C be n rows (one for each neuron) and two columns (one for the response to the latent input)
for mm = 1:yDim
    C(mm,1) = 2 + 0.00541*exp((0.049)*mm);      %Set the linear equation for how C changes for each neuron for input conductance
    C(mm,2) = 0.5 + 0.001666 * exp((0.048)*mm); %Set linear equaiton for C for AHP duration
end


LinMult = AboveRTE ./5;                          %Create column array of force differences to model sigmoid relation between input force and discharge rate. 0.5 is just a scaling factor
LinMult = tanh(LinMult) * 0.5 + 0.5;             %Set up an approximate sigmoid relation
LinMult = LinMult(1:ActiveUnits(500),:);         %Get to correct size for a given trial

FiringRate = zeros(ActiveUnits(500),500);           %preallocate a discharge rate (intensity) matrix
AboveRTE = AboveRTE(1:ActiveUnits(500),1:500);
NoiseMat = ones(ActiveUnits(500),1); 
NoiseMat = NoiseMat .* (AboveRTE(:,500).* 5);       %In this way, a unit closer to recruitment threshold receives more noise than one further away
NoiseMat(NoiseMat<0) = 0;                           %5 is another scaling factor


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

figure(1);
plotSpikeRaster(logical(dsSub),'PlotType','Scatter');
set(gcf,'color','white')



plot(params.model.Xpca(1,:),'k', 'LineWidth',2)
hold on
%plot(params.model.Xpca(2,:), 'LineWidth',2, 'Color',[0.501960813999176 0.501960813999176 0.501960813999176])
set(gca,'TickDir','out', 'fontsize',16);%, 'XTickLabel',{'0','1','2','3','4','5'}); 
box off
set(gcf,'color','w')
hold off
set(gcf, 'WindowStyle','normal')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 4], 'PaperUnits', 'Inches', 'PaperSize', [4, 6])