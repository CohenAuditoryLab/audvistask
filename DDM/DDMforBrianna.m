function DDMforBrianna
% function OutDT=DDMforBrianna

% 1st column: percentage of high tone (I mapped this percentage to signed coherence like this: (proportion-50).*2)
% 2nd column: monkey's choice (1: high-tone choice, 0: low-tone choice)
% 3rd column: response time (sec)

close all

disp('change here to specify the directory for the data')
%cd ('C:\JT\ProcessedData\JoshPFC');
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/DDM/');
pause


load PopBehavior
% This is behavioral data from monkey Tyr in PFC recroding sessions.
% PFC sessions
% 'tyr120524.mat';  
% 'tyr120529.mat';  
% 'tyr120531.mat';  
% 'tyr120604.mat';  
% 'tyr120607.mat';  
% 'tyr120612.mat';  
% 'tyr120614.mat';  
% 'tyr120615.mat';  
% 'tyr120626.mat';  
% 'tyr120627.mat';  
% 'tyr130909.mat';  
% 'tyr130917.mat';  
% 'tyr130926.mat';  
% 'tyr130927.mat';  
% 'tyr131010.mat';  
% 'tyr131224.mat';  
% 'tyr131226.mat';  
% 'tyr140106.mat';  
% 'tyr140117.mat';  
% 'tyr140128.mat';  
% 'tyr140131.mat';  
% 'tyr140204.mat';  
% 'tyr140213.mat';  
% 'tyr140227.mat';  
% 'tyr140228.mat';  
% 'tyr140325.mat';  
% 'tyr140327.mat';  
% 'tyr140331.mat';  
% 'tyr140407.mat';

%coherence bins
cbins = [ ...
    -100  -99
    -99  -60
    -60  -20
    -20    20
    20   60
    60  100
    100  101];

nbins = size(cbins,1);
cax   = mean(cbins,2);
cfax  = -100:.1:100;
pmf   = NaN(nbins,1);
cmf   = NaN(nbins,2);

% X0  = [200 200 200 200 200];
% Xlb = [0.01 0.01 0.01 0.01 0.01];
% Xub = [50000 50000 50000 50000 50000];

Behavior=PopBehavior;
disp('change here to specify the directory for matlab codes')
%cd ('C:\JT\MATLAB\Josh\DDM3');
cd ('/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/DDM/');
pause

ntrials = size(Behavior,1);
Lch     = Behavior(:,2)==1;

% data matrix, columns are:
%   1. re-scale signal strength to [-100, 100]
%   2. choice (1: high-tone choice, 0: low-tone choice)
%   3. RT (sec)
%   4. correct (1) or error (0)
scoh = Behavior(:,1)*2-100;
Lcor = (Behavior(:,1)>= -20 & Behavior(:,1)<20) | (Behavior(:,1)>50&Lch) | (Behavior(:,1)<50&~Lch);
data = cat(2, scoh, Behavior(:,2), Behavior(:,3).*1000, double(Lcor));

% make selection array, compute pmf cmf
Lcoh  = false(ntrials, nbins);
for cc = 1:nbins
    Lcoh(:,cc) = scoh>=cbins(cc,1) & scoh<cbins(cc,2);
    cax(cc) = mean(scoh(Lcoh(:,cc)));
    
    pmf(cc) = sum(Behavior(Lcoh(:,cc)&Lch,2))./sum(Lcoh(:,cc)).*100;
    
    cmf(cc,1) = nanmean(data(Lcoh(:,cc)& Lch,3));
    cmf(cc,2) = nanmean(data(Lcoh(:,cc)&~Lch,3));
end

% RT for -100% decision B
RThigh=cmf(end,1);%decisionA
RTlow=cmf(1,2);%decisionB
X0  = [200 200 200 200 200];
Xlb = [0.01 0.01 0.01 0.01 0.01];
%     Xub = [50000 50000 50000 RThigh-170 RTlow-170];
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

% save the best
if errg < err0
    fits = fitsg;
end

fdat(1,:) = [fits lapse];

%     %% PLOTZ
[ps,rts] = fitJT_val_simple5L(cfax, fits, lapse);

subplot(3,2,1); cla reset; hold on;
plot(cax, pmf, 'k.');
plot(cfax, ps.*100, 'r-');
title(sprintf('%.2f, %.2f, %.2f, %.2f, %.2f', ...
    fits(1), fits(2), fits(3), fits(4), fits(5)))
plot([-100 0], [lapse*100 lapse*100], 'k--');
plot([0 100], [100-lapse*100 100-lapse*100], 'k--');
xlabel('Coherence (%): +100 means all high tones')
ylabel('high-tone choice (%)')

subplot(3,2,3); cla reset; hold on;
%     plot(cax(cax>=0), cmf(cax>=0,1), 'k.');
%     plot(cax(cax<=0), cmf(cax<=0,2), 'b.');
plot(cax(cax>=-20), cmf(cax>=-20,1), 'k.');
plot(cax(cax<=20), cmf(cax<=20,2), 'b.');
plot(cfax, rts, 'r-');
plot([0 100], fits([4 4]), 'k--');
plot([-100 0], fits([5 5]), 'b--');
xlabel('Coherence (%): +100 means all high tones')
ylabel('Response time (ms)')

subplot(3,2,5); cla reset; hold on;
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


