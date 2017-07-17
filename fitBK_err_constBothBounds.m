function err_ = fitBK_err_constBothBounds(fits, data, lapse)
% 
% fits are from fitBK_val_constBounds3L

l = size(data, 2);

cohs = [data(:,1), data(:,5), data(:,9)];

if l >= 13
    cohs = [cohs, data(:, 13)];
end 
if l >= 17
    cohs = [cohs, data(:, 17)];
end

[ps, rts] = fitBK_val_constBounds3L(cohs, fits, lapse);

% logL of pmf
logpPMF1 = log(binopdf(data(:,2), 1, ps(:,1)));
logpPMF1(~isfinite(logpPMF1)) = -200;
logLp1   = sum(logpPMF1);

logpPMF2 = log(binopdf(data(:,6), 1, ps(:, 2)));
logpPMF2(~isfinite(logpPMF2)) = -200;
logLp2   = sum(logpPMF2);

logpPMF3 = log(binopdf(data(:,10), 1, ps(:, 3)));
logpPMF3(~isfinite(logpPMF3)) = -200;
logLp3   = sum(logpPMF3);

if l >= 13
    logpPMF4 = log(binopdf(data(:,14), 1, ps(:, 4)));
    logpPMF4(~isfinite(logpPMF4)) = -200;
    logLp4   = sum(logpPMF4);
end

if l >= 17
    logpPMF5 = log(binopdf(data(:,18), 1, ps(:, 5)));
    logpPMF5(~isfinite(logpPMF5)) = -200;
    logLp5   = sum(logpPMF5);
end 

% logL of cmf -- only for non-zero CORRECT trials
% assume RT variance is the same per condition, measured from the data
Lgood1   = data(:,4) == 1;
rts1 = rts(:,1);
msrts1   = data(Lgood1,3) - rts1(Lgood1);
logpCMF1 = log(normpdf(msrts1, 0, std(msrts1)));
logpCMF1(~isfinite(logpCMF1)) = -200; 
logLc1   = sum(logpCMF1);

Lgood2   = data(:,8) == 1;
rts2 = rts(:, 2);
msrts2   = data(Lgood2,3) - rts2(Lgood2);
logpCMF2 = log(normpdf(msrts2, 0, std(msrts2)));
logpCMF2(~isfinite(logpCMF2)) = -200; 
logLc2   = sum(logpCMF2);

Lgood3   = data(:,12) == 1;
rts3 = rts(:, 3);
msrts3   = data(Lgood3,3) - rts3(Lgood3);
logpCMF3 = log(normpdf(msrts3, 0, std(msrts3)));
logpCMF3(~isfinite(logpCMF3)) = -200; 
logLc3   = sum(logpCMF3);

if l >= 13
    Lgood4   = data(:,16) == 1;
    rts4 = rts(:, 4);
    msrts4   = data(Lgood4,3) - rts4(Lgood4);
    logpCMF4 = log(normpdf(msrts4, 0, std(msrts4)));
    logpCMF4(~isfinite(logpCMF4)) = -200;
    logLc4   = sum(logpCMF4);
end 

if l >= 17
    Lgood5   = data(:,20) == 1;
    rts5 = rts(:, 5);
    msrts5   = data(Lgood5,3) - rts5(Lgood5);
    logpCMF5 = log(normpdf(msrts5, 0, std(msrts5)));
    logpCMF5(~isfinite(logpCMF5)) = -200;
    logLc5   = sum(logpCMF5);
end

% err_ is negative sum of the log likelihoods
err_ = -(logLp1 + logLc1 + logLp2 + logLc2 + logLp3 + logLc3);
if l > 12
    err_ = err_ -(logLp4 + logLc4);
end
if l >= 17
    err_ = err_ -(logLp5 + logLc5);
end
if isnan(err_)
    err_ = inf;
end

end

