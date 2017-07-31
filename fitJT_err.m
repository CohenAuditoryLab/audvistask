function err_ = fitJT_err(fits, data, lapse)
% 
% fits are from fitJT_val_simple5L
% 
% data columns are:
% fits are val* function specific
% ivar is matrix of independent variables
% data is matrix, rows correspond to independent variables, columns are:
%   1 .. # T1 choices
%   2 .. # correct
%   3 .. # total
%   4 .. mean rt correct
%   5 .. var rt correct
%   6 .. mean rt error
%   7 .. var rt error
% Generic error function that uses handles to compute
%   the predictions and the error

% Each file contains data for each session, matrix name = Behavior 
% 1st column: signed coh
% 2nd column: monkey's choice (1: high-tone choice, 0: low-tone choice) 
% 3rd column: response time (sec)
% 4th column: correct

[ps, rts] = fitJT_val_simple5L(data(:,1), fits, lapse);

% logL of pmf
logpPMF = log(binopdf(data(:,2), 1, ps));
logpPMF(~isfinite(logpPMF)) = -200;
logLp   = sum(logpPMF);

% logL of cmf -- only for non-zero CORRECT trials
% assume RT variance is the same per condition, measured from the data
Lgood   = data(:,4) == 1;
msrts   = data(Lgood,3) - rts(Lgood);
logpCMF = log(normpdf(msrts, 0, std(msrts)));
logpCMF(~isfinite(logpCMF)) = -200; % if it's not Scottish it's CRAP
logLc   = sum(logpCMF);

% err_ is negative sum of the log likelihoods
err_ = -(logLp + logLc);
if isnan(err_)
    err_ = inf;
end
