function MyDDMAuditory
% Data File being used:
% 1st column: percentage of high tone (mapped this percentage to signed coherence like this: (proportion-50).*2)
% 2nd column: monkey's choice (1: high-tone choice, 0: low-tone choice)
% 3rd column: response time (sec)
% 4th column: success of choice (1: correct, 0: incorrect)

close all

%cd into file that holds data 
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/');
%load in data
PopBehavior = csvread('DDM_AudVisTask_v1_Brianna_170530_1345.csv', 1, 0);
%csvread('SampleData.csv', 0, 0); 

%coherence bins
cbins = [ ...
    -100  -99
    -99  -60
    -60  -20
    -20    20
    20   60
    60  100
    100  101];
% cbins = [...
%         -100 -75
%         -75 -50
%         -50 -25 
%         -25 25
%         25 50
%         50 75
%         75 100];

%calculate number of coherence bins 
nbins = size(cbins,1);
%calculate average value of each coherence bin - output vector 
cax   = mean(cbins,2);
%create vector from -100 to 100
cfax  = -100:.1:100;
%initialize vectors to hold pmf and cmf values
pmf   = NaN(nbins,1);
cmf   = NaN(nbins,2);

%set data to variable
Behavior=PopBehavior;
%enter file for fit functions
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/DDM/');

%calculate number of trials 
ntrials = size(Behavior,1);
%create vector of true = high false = low; ensures all values are 0 or 1 
Lch     = Behavior(:,2)==1;

% create a data matrix, columns are:
%   1. re-scale signal strength to [-100, 100]
%   2. choice (1: high-tone choice, 0: low-tone choice)
%   3. RT (ms)
%   4. correct (1) or error (0)
scoh = Behavior(:,1)*200-100;
Lcor = Behavior(:,4); %(Behavior(:,1)>= -20 & Behavior(:,1)<20) | (Behavior(:,1)>50&Lch) | (Behavior(:,1)<50&~Lch);
data = cat(2, scoh, Behavior(:,2), Behavior(:,3), double(Lcor));

% make selection array, compute pmf and cmf
% PMF - a function that gives the probability that a discrete random variable is exactly equal to some value
% CMF -  the probability that X will take a value less than or equal to x
Lcoh  = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh(:,cc) = scoh>=cbins(cc,1) & scoh<cbins(cc,2);
    cax(cc) = mean(scoh(Lcoh(:,cc)));
    
    pmf(cc) = sum(Behavior(Lcoh(:,cc)&Lch,2))./sum(Lcoh(:,cc)).*100;
    
    cmf(cc,1) = nanmean(data(Lcoh(:,cc)& Lch,3));
    cmf(cc,2) = nanmean(data(Lcoh(:,cc)&~Lch,3));
end

% RT for -100% decision B and 100% decision A
RThigh = cmf(end,1); % decision A
RTlow = cmf(1,2); % decision B
X0  = [200 200 200 200 200];
Xlb = [0.01 0.01 0.01 0.01 0.01];
%Xub = [50000 50000 50000 RThigh-170 RTlow-170];
Xub = [50000 50000 50000 2000 2000];

% compute lapse
Llapse = abs(scoh)>=90;
lapse = 1-sum(Llapse&Lcor)./sum(Llapse);

% get err from initial fit
err0 = fitJT_err(X0, data, lapse);

% fit it using pattern search
[fits,err] = patternsearch(@(x)fitJT_err(x, data, lapse), ...
    X0, [], [], [], [], Xlb, Xub, [], ...
    psoptimset('MaxIter', 5000, 'MaxFunEvals', 5000));

% possibly seed as init values
if err < err0
    X0g  = fits;
    err0 = err;
else
    X0g = X0;
end

% now gradient descent
[fitsg,errg] = fmincon(@(x)fitJT_err(x, data, lapse), ...
    X0g, [], [], [], [], Xlb, Xub, [], ...
    optimset('Algorithm', 'active-set', ...
    'MaxIter', 30000, 'MaxFunEvals', 30000));%, 'Display', 'iter'));

% save the best between the two methods
if errg < err0
    fits = fitsg;
end

fdat(1,:) = [fits lapse];

%%% PLOTZ
[ps,rts] = fitJT_val_simple5L(cfax, fits, lapse);

subplot(3,1,1); cla reset; hold on;
plot(cax, pmf, 'k.', 'MarkerSize', 8);
plot(cfax, ps.*100, 'r-');
title(sprintf('%.2f, %.2f, %.2f, %.2f, %.2f', ...
    fits(1), fits(2), fits(3), fits(4), fits(5)))
plot([-100 0], [lapse*100 lapse*100], 'b--');
plot([0 100], [100-lapse*100 100-lapse*100], 'k--');
ylim([-1 101])
xlabel('Coherence (%): +100 means all high tones')
ylabel('high-tone choice (%)')

subplot(3,1,2); cla reset; hold on;
plot(cax(cax>=0), cmf(cax>=0,1), 'k.');
plot(cax(cax<=0), cmf(cax<=0,2), 'b.');
%plot(cax(cax>=-20), cmf(cax>=-20,1), 'k.', 'MarkerSize', 8);
%plot(cax(cax<=20), cmf(cax<=20,2), 'b.', 'MarkerSize', 8);
plot(cfax, rts, 'r-');
plot([0 100], fits([4 4]), 'k--');
plot([-100 0], fits([5 5]), 'b--');
ylim([-50 2000])
xlabel('Coherence (%): +100 means all high tones')
ylabel('Response time (ms)')

subplot(3,1,3); cla reset; hold on;
plot(cfax(cfax<0), rts(cfax<0)-fits(5), 'r-');hold on;
plot(cfax(cfax>0), rts(cfax>0)-fits(4), 'r-');
xlabel('Coherence (%): +100 means all high tones')
ylabel('Decision time (ms): RT-nonDT')

%Decision time
% OutDT=[rts(cfax<0)-fits(5) min([max(rts(cfax < 0)-fits(5)) max(rts(cfax > 0)-fits(4))]) rts(cfax>0)-fits(4)];

%% Logistic fit
% if you want to test with logistic fit, run below.
% PsyData=[ones(size(PopBehavior,1),1) PopBehavior(:,1:2)];
% [fits_,sems_,stats_,preds_,resids_] = logist_fit(PsyData, 'lug',[],lapse);
% X=[ones(2000,1) linspace(0, 100, 2000)'];
% ys = 100*feval(@logist_valLUg, fits_,X,lapse);
% subplot(3,1,1); 
% plot(linspace(-100, 100, 2000)',ys,'g');


