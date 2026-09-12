function [output] = NGAR(data, target, mumean1b, mulambdastar, burnin, numbofits, every, n_condition, p_condition, options)

%NGAR Dynamic regression with optional normal/Poisson count approximations.
%   OUT = NGAR(X,Y,MUMEAN1B,MULAMBDASTAR,BURNIN,NUMBOFITS,EVERY) uses
%   MATLAB count samplers.
%   OUT = NGAR(...,N_CONDITION,P_CONDITION) uses the approximation rules.
%   They provide an alternative when ordinary large-count draws become
%   computationally impractical. Set N_CONDITION=Inf to disable them.
%   P_CONDITION is the half-width around 0.5 for the normal-binomial branch.
%   X is T-by-p; Y may be a row or column vector. Include an intercept in X
%   explicitly: the first column receives the original special prior.
%
%   The approximation option uses approximate proposal draws; its invariant
%   distribution has not been established as identical to the ordinary route.
%   The NGAR parameter bounds, adaptation schedule and
%   sampling calculations are retained. See docs/DISPLAY_DIAGNOSTICS.md.
%
%   Adapted from the MATLAB code accompanying Kalli and Griffin (2014),
%   Journal of Econometrics 178, 779-793,
%   doi:10.1016/j.jeconom.2013.10.012. Carl Gower added the optional count
%   approximations. Published on GitHub with written permission from Maria
%   Kalli and Jim E. Griffin. No named standard software licence is asserted;
%   see PROVENANCE.md and LICENSE.md.
%   Requires Statistics and Machine Learning Toolbox and kf_NGAR.m.
%   Set rng before calling when reproducibility is required.
%
%   OUT = NGAR(...,N_CONDITION,P_CONDITION,OPTIONS) accepts a scalar struct:
%     OPTIONS.PrintEvery     Completed iterations between displays (default 20).
%                            Set to 0 for no progress display.
%     OPTIONS.SaveParameters Parameters whose draws are saved (default 'all').
%                            Example: {'beta','Psi'}. Use {} for diagnostics only.
%   With ordinary count samplers, OPTIONS may instead be the eighth input.
%   All parameters are sampled; SaveParameters controls storage and output only.
%   OUT.diagnostics records settings, timing, RNG states and count-branch usage.

narginchk(7, 10);
if nargin < 10
    options = struct;
end
if nargin == 8
    if ~isstruct(n_condition)
        error('NGAR:MissingThreshold', 'Supply both thresholds, or an options struct as the eighth input.');
    end
    options = n_condition;
    n_condition = Inf;
    p_condition = 0.25;
elseif nargin == 7
    n_condition = Inf;
    p_condition = 0.25;
end
validateattributes(data, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, mfilename, 'data');
[T, p] = size(data);
if T < 2
    error('NGAR:TooFewObservations', 'At least two observations are required.');
end
validateattributes(target, {'numeric'}, {'vector', 'real', 'finite', 'numel', T}, mfilename, 'target');
validateattributes(mumean1b, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'mumean1b');
validateattributes(mulambdastar, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, mfilename, 'mulambdastar');
validateattributes(burnin, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, mfilename, 'burnin');
validateattributes(numbofits, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'numbofits');
validateattributes(every, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'every');
validateattributes(n_condition, {'numeric'}, {'scalar', 'real', 'nonnan', '>=', 1}, mfilename, 'n_condition');
validateattributes(p_condition, {'numeric'}, {'scalar', 'real', 'finite', '>=', 0, '<=', 0.5}, mfilename, 'p_condition');
validateattributes(options, {'struct'}, {'scalar'}, mfilename, 'options');
unknownOptions = setdiff(fieldnames(options), {'PrintEvery'; 'SaveParameters'});
if ~isempty(unknownOptions)
    error('NGAR:UnknownOption', 'Unknown option: %s.', unknownOptions{1});
end
if ~isfield(options, 'PrintEvery')
    options.PrintEvery = 20;
end
if ~isfield(options, 'SaveParameters')
    options.SaveParameters = 'all';
end
validateattributes(options.PrintEvery, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, mfilename, 'PrintEvery');
parameterNames = {'beta', 'sigmasq', 'Psi', 'lambda', 'mu', 'rhobeta', 'rho', ...
    'lambdasigma', 'musigma', 'rhosigma', 'lambdastar', 'mustar'};
requestedParameters = options.SaveParameters;
if ischar(requestedParameters) && isrow(requestedParameters)
    requestedParameters = {requestedParameters};
elseif isstring(requestedParameters) && (isvector(requestedParameters) || isempty(requestedParameters))
    requestedParameters = cellstr(requestedParameters(:));
end
if ~iscellstr(requestedParameters) || ~(isvector(requestedParameters) || isempty(requestedParameters))
    error('NGAR:InvalidSaveParameters', 'SaveParameters must be ''all'' or a list of parameter names.');
end
requestedParameters = strtrim(requestedParameters(:)');
if isscalar(requestedParameters) && strcmpi(requestedParameters{1}, 'all')
    requestedParameters = parameterNames;
end
saveParameters = struct;
for parameterIndex = 1:numel(parameterNames)
    parameterName = parameterNames{parameterIndex};
    saveParameters.(parameterName) = any(strcmpi(parameterName, requestedParameters));
end
for parameterIndex = 1:numel(requestedParameters)
    if ~any(strcmpi(requestedParameters{parameterIndex}, parameterNames))
        error('NGAR:InvalidSaveParameters', 'Unknown parameter: %s.', requestedParameters{parameterIndex});
    end
end
options.SaveParameters = parameterNames(cellfun(@(name) saveParameters.(name), parameterNames));
inputTargetSize = size(target);
runTimer = tic;
initialRng = rng;
countUsage = zeros(1, 7);
invalidPathProposals = 0;

data = double(data);
target = double(target(:)'); % Preserve the row-vector convention in the original residual calculation.

np_threshold = 10;
limit1 = 0.9999;
lambda_min = 0.1;
lambda_max = 100;
mu_min = 0.01;
mu_max = 100;

rho = 0.97 * ones(1, p);
rhobeta = 0.97 * ones(1, p);
lambda = min(lambda_max, max(lambda_min, mulambdastar)) * ones(1, p);
mu = min(mu_max, max(mu_min, mumean1b)) * ones(1, p);
Delta = rho ./ (1 - rho) .* lambda ./ mu;

start_samples = 500;
start_adap = 1000;

newbeta = zeros(T, p);
Psi = ones(T, p);
kappa = zeros(T, p);
for j = 1:p
    kappa(:, j) = poissrnd(Delta(j) * Psi(:, j));
    countUsage(5) = countUsage(5) + T;
end
kappasigmasq = ones(T, 1);
lambdasigma = 3;
musigma = 0.03;
rhosigma = 0.95;

sigmasq = musigma * ones(T, 1);

lambdastar = 1;
mustar = 1;

Psisd = 0.01 * ones(T, p);
loglambdasd = log(0.1) * ones(1, p);
logmeansd = log(0.1) * ones(1, p);
logscale = 0.5 * log(2.4^2 / 4) * ones(1, p);
logrhobetasd = log(0.1) * ones(1, p);
logrhosd = log(0.1) * ones(1, p);
logrhosigmasd = log(0.01);
loglambdasigmasd = log(0.001);
loggammasigmasqsd = log(0.001);
logscalesigmasq = 0.5 * log(2.4^2 / 3);
mugammasd = 0.003;
vgammasd = 0.003;
logsigmasqsd = log(0.001) * ones(1, T);
logkappaq = 4 / 3 * ones(T - 1, p);
logkappasigmasqq = 4 / 3 * ones(T, 1);

kappaaccept = 0;
kappacount = 0;
kappalambdasigmaccept = 0;
kappasigmasqcount = 0;
sigmasqparamaccept = 0;
sigmasqparamcount = 0;
Psiparamaccept = zeros(1, p);
Psiparamcount = zeros(1, p);
mugammaaccept = 0;
mugammacount = 0;
vgammaaccept = 0;
vgammacount = 0;
sigmasq1accept = zeros(1, T);
sigmasq1count = zeros(1, T);
numberofiterations = burnin + every * numbofits;

% Allocate histories only for requested parameters (otherwise the last dimension is zero).
holdPsi = zeros(T, p, numbofits * saveParameters.Psi);
holdbeta = zeros(T, p, numbofits * saveParameters.beta);
holdsigmasq = zeros(T, numbofits * saveParameters.sigmasq);
holdlambda = zeros(p, numbofits * saveParameters.lambda);
holdmu = zeros(p, numbofits * saveParameters.mu);
holdrhobeta = zeros(p, numbofits * saveParameters.rhobeta);
holdrho = zeros(p, numbofits * saveParameters.rho);
holdlambdasigma = zeros(1, numbofits * saveParameters.lambdasigma);
holdmusigma = zeros(1, numbofits * saveParameters.musigma);
holdrhosigma = zeros(1, numbofits * saveParameters.rhosigma);
holdlambdastar = zeros(1, numbofits * saveParameters.lambdastar);
holdmustar = zeros(1, numbofits * saveParameters.mustar);

sum1 = zeros(4, p);
sum2 = zeros(10, p);
sum1sigmasq = zeros(3, 1);
sum2sigmasq = zeros(6, 1);

[meankf, varkf, ~] = kf_NGAR(data, target, Psi, sigmasq, rhobeta);

cholstar = chol(varkf(:, :, T))';
newbeta(T, :) = (meankf(:, T) + cholstar * randn(size(cholstar, 2), 1))';
for i = (T-1):-1:1
    Gkal = diag(rhobeta .* sqrt(Psi(i + 1, :) ./ Psi(i, :)));
    invQ = diag(1 ./ (1 - rhobeta.^2) .* 1./Psi(i + 1, :));
    % Scale the matrix to unit diagonal, invert it, then undo the scaling.
    covScale = sqrt(diag(varkf(:, :, i)));
    invvarkf = ((varkf(:, :, i) ./ (covScale * covScale')) \ eye(numel(covScale))) ./ (covScale * covScale');
    precisionfb = invvarkf + Gkal' * invQ * Gkal;
    precisionScale = sqrt(diag(precisionfb));
    varfb = ((precisionfb ./ (precisionScale * precisionScale')) \ eye(numel(precisionScale))) ./ (precisionScale * precisionScale');
    meanfb = varfb * (invvarkf * meankf(:, i) + Gkal' * invQ * newbeta(i + 1, :)');
    cholstar = chol(varfb)';
    newbeta(i, :) = (meanfb + cholstar * randn(size(cholstar, 2), 1))';
end

beta = newbeta(:, 1:p);

checkstar = 1;

for it = 1:numberofiterations

    if ( checkstar == 1 )

        % Update coefficient scales.
        for i = 1:T
            for j = 1:p
                newPsi = Psi(i, j) * exp(Psisd(i, j) * randn);

                if (i == 1 )
                    loglike = (lambda(j) - 1) * log(Psi(1, j)) - lambda(j) * Psi(1, j) / mu(j);
                    pnmean = Psi(i, j) * Delta(j);
                    loglike = loglike - pnmean + kappa(i, j) * log(pnmean);
                    loglike = loglike - 0.5 * log(Psi(1, j)) - 0.5 * beta(1, j)^2 / Psi(1, j);
                    var1 = Psi(2, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(2, j) / Psi(1, j)) * beta(1, j);
                    loglike = loglike - 0.5 * (beta(2, j) - mean1)^2 / var1;

                    newloglike = (lambda(j) - 1) * log(newPsi) - lambda(j) * newPsi / mu(j);
                    pnmean = newPsi * Delta(j);
                    newloglike = newloglike - pnmean + kappa(i, j) * log(pnmean);
                    newloglike = newloglike - 0.5 * log(newPsi) - 0.5 * beta(1, j)^2 / newPsi;
                    var1 = Psi(2, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(2, j) / newPsi) * beta(1, j);
                    newloglike = newloglike - 0.5 * (beta(2, j) - mean1)^2 / var1;
                elseif ( i < T )
                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    loglike = (lam1 - 1) * log(Psi(i, j)) - gam1 * Psi(i, j);
                    pnmean = Psi(i, j) * Delta(j);
                    loglike = loglike - pnmean + kappa(i, j) * log(pnmean);
                    var1 = Psi(i, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i, j) / Psi(i - 1, j)) * beta(i - 1, j);
                    loglike = loglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                    var1 = Psi(i + 1, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i + 1, j) / Psi(i, j)) * beta(i, j);
                    loglike = loglike - 0.5 * (beta(i + 1, j) - mean1)^2 / var1;

                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    newloglike = (lam1 - 1) * log(newPsi) - gam1 * newPsi;
                    pnmean = newPsi * Delta(j);
                    newloglike = newloglike - pnmean + kappa(i, j) * log(pnmean);
                    var1 = newPsi * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(newPsi / Psi(i - 1, j)) * beta(i - 1, j);
                    newloglike = newloglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                    var1 = Psi(i + 1, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i + 1, j) / newPsi) * beta(i, j);
                    newloglike = newloglike - 0.5 * (beta(i + 1, j) - mean1)^2 / var1;
                else
                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    loglike = (lam1 - 1) * log(Psi(i, j)) - gam1 * Psi(i, j);
                    var1 = Psi(i, j) * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(Psi(i, j) / Psi(i - 1, j)) * beta(i - 1, j);
                    loglike = loglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;

                    lam1 = lambda(j) + kappa(i - 1, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    newloglike = (lam1 - 1) * log(newPsi) - gam1 * newPsi;
                    var1 = newPsi * (1 - rhobeta(j)^2);
                    mean1 = rhobeta(j) * sqrt(newPsi / Psi(i - 1, j)) * beta(i - 1, j);
                    newloglike = newloglike - 0.5 * log(var1) - 0.5 * (beta(i, j) - mean1)^2 / var1;
                end


                logaccept = newloglike - loglike + log(newPsi) - log(Psi(i, j));
                accept = 1;
                if ( isnan(logaccept) || isinf(logaccept) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end

                Psisd(i, j) = Psisd(i, j) + (accept - 0.3) / it^0.6;

                if ( rand < accept )
                    Psi(i, j) = newPsi;
                end
            end
        end



        for i = 1:(T-1)
            for j = 1:p
                newkappa = kappa(i, j) + (2 * (rand < 0.5) - 1) * geornd(1 / (1 + exp(logkappaq(i, j))));

                if ( newkappa < 0 )
                    accept = 0;
                else
                    lam1 = lambda(j) + kappa(i, j);
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    loglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(Psi(i + 1, j));
                    pnmean = Psi(i, j) * Delta(j);
                    loglike = loglike + kappa(i, j) * log(pnmean) - gammaln(kappa(i, j) + 1);

                    lam1 = lambda(j) + newkappa;
                    gam1 = lambda(j) / mu(j) + Delta(j);
                    newloglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(Psi(i + 1, j));
                    pnmean = Psi(i, j) * Delta(j);
                    newloglike = newloglike + newkappa * log(pnmean) - gammaln(newkappa + 1);

                    logaccept = newloglike - loglike;

                    accept = 1;
                    if ( isnan(logaccept) || isinf(logaccept) )
                        accept = 0;
                    elseif ( logaccept < 0 )
                        accept = exp(logaccept);
                    end
                end

                kappaaccept = kappaaccept + accept;
                kappacount = kappacount + 1;

                if ( rand < accept )
                    kappa(i, j) = newkappa;
                end

                logkappaq(i, j) = logkappaq(i, j) + 1 / it^0.55 * (accept - 0.3);

                if ( ~isfinite(kappa(i, j)) || ~isreal(kappa(i, j)) )
                    error('NGAR:InvalidCount', 'An auxiliary count is nonfinite or complex.');
                end
            end
        end


        for i = 1:T
            chi1 = (target(i) - sum(data(i, :) .* beta(i, :)))^2;
            if ( i == 1 )
                lam1 = kappasigmasq(1) + lambdasigma - 0.5;
                psi1 = 2 * (lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma);
            elseif ( i == T )
                lam1 = kappasigmasq(i - 1) + lambdasigma - 0.5;
                psi1 = 2 * (lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma);
            else
                lam1 = kappasigmasq(i) + kappasigmasq(i - 1) + lambdasigma - 0.5;
                psi1 = 2 * (lambdasigma / musigma + 2 * rhosigma / (1 - rhosigma) * lambdasigma / musigma);
            end

            newsigmasq = sigmasq(i) * exp(exp(logsigmasqsd(i)) * randn);

            loglike = (lam1 - 1) * log(sigmasq(i)) - 0.5 * chi1 ./ sigmasq(i) - 0.5 * psi1 * sigmasq(i);
            newloglike = (lam1 - 1) * log(newsigmasq) - 0.5 * chi1 ./ newsigmasq - 0.5 * psi1 * newsigmasq;

            logaccept = newloglike - loglike + log(newsigmasq) - log(sigmasq(i));

            accept = 1;
            if ( isnan(logaccept) || isinf(logaccept) )
                accept = 0;
            elseif ( logaccept < 0 )
                accept = exp(logaccept);
            end

            sigmasq1accept(i) = sigmasq1accept(i) + accept;
            sigmasq1count(i) = sigmasq1count(i) + 1;

            if ( rand < accept )
                sigmasq(i) = newsigmasq;
            end

            logsigmasqsd(i) = logsigmasqsd(i) + 1 / it^0.55 * (accept - 0.3);
        end

        for i = 1:(T-1)
            newkappasigmasq = kappasigmasq(i) + (2 * (rand < 0.5) - 1) * geornd(1 / (1 + exp(logkappasigmasqq(i))));

            if ( newkappasigmasq < 0 )
                accept = 0;
            else
                lam1 = lambdasigma + kappasigmasq(i);
                gam1 = lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                loglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(sigmasq(i + 1));
                pnmean = sigmasq(i) * rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                loglike = loglike + kappasigmasq(i) * log(pnmean) - gammaln(kappasigmasq(i) + 1);

                lam1 = lambdasigma + newkappasigmasq;
                gam1 = lambdasigma / musigma + rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                newloglike = lam1 * log(gam1) - gammaln(lam1) + (lam1 - 1) * log(sigmasq(i + 1));
                pnmean = sigmasq(i) * rhosigma / (1 - rhosigma) * lambdasigma / musigma;
                newloglike = newloglike + newkappasigmasq * log(pnmean) - gammaln(newkappasigmasq + 1);

                logaccept = newloglike - loglike;

                accept = 1;
                if ( isnan(logaccept) || isinf(logaccept) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end
            end

            kappalambdasigmaccept = kappalambdasigmaccept + accept;
            kappasigmasqcount = kappasigmasqcount + 1;

            if ( rand < accept )
                kappasigmasq(i) = newkappasigmasq;
            end

            logkappasigmasqq(i) = logkappasigmasqq(i) + 1 / it^0.55 * (accept - 0.3);


            if ( ~isfinite(kappasigmasq(i)) || ~isreal(kappasigmasq(i)) )
                error('NGAR:InvalidCount', 'An auxiliary count is nonfinite or complex.');
            end
        end
    end

    % Select the coefficient block for the joint parameter/path updates.
    zstar = rand(1, p) < (5 / p);


    targetstar = target - (sum(data(:, zstar==0) .* beta(:, zstar==0), 2))';
    datastar = data(:, zstar==1);
    Psistar = Psi(:, zstar==1);
    kappastar = kappa(:, zstar==1);

    [meankf, varkf, loglike] = kf_NGAR(datastar, targetstar, Psistar, sigmasq, rhobeta(zstar==1));




    xstar = [log(lambdasigma); log(musigma); log(rhosigma) - log(1 - rhosigma)];

    if ( it < 100 )
        newxstar = xstar + [exp(loglambdasigmasd); exp(loggammasigmasqsd); exp(logrhosigmasd)] .* randn(3, 1);
    else
        varstar1 = ([sum2sigmasq(1) sum2sigmasq(2) sum2sigmasq(4); sum2sigmasq(2) sum2sigmasq(3) sum2sigmasq(5); sum2sigmasq(4) sum2sigmasq(5) sum2sigmasq(6)] - sum1sigmasq * sum1sigmasq' / (it - 1)) / (it - 2);

        [proposalChol, proposalStatus] = chol(exp(logscalesigmasq) * varstar1);
        if proposalStatus == 0
            newxstar = xstar + proposalChol' * randn(3, 1);
        else
            % If the learned covariance cannot be factored, use independent parameter steps.
            newxstar = xstar + [exp(loglambdasigmasd); exp(loggammasigmasqsd); exp(logrhosigmasd)] .* randn(3, 1);
        end
    end

    newlambdasigma = exp(newxstar(1));
    newmusigma = exp(newxstar(2));
    newrhosigma = exp(newxstar(3)) / (1 + exp(newxstar(3)));

    if ( ~all(isfinite([newlambdasigma, newmusigma, newrhosigma])) || newrhosigma <= 0 || newrhosigma > limit1 || ...
        (newlambdasigma < lambda_min) || (newlambdasigma > lambda_max ) || ...
        (newmusigma         < mu_min    ) || (newmusigma     > mu_max     ))
        accept = 0;
    else
        validPath = true;
        newsigmasq = sigmasq;
        newkappasigmasq = kappasigmasq;
        newsigmasq(1) = sigmasq(1) * (lambdasigma / newlambdasigma) * newmusigma / musigma;

        if ( newlambdasigma > lambdasigma )
            newsigmasq(1) = newsigmasq(1) + gamrnd(newlambdasigma - lambdasigma, newmusigma / newlambdasigma);
        elseif newlambdasigma < lambdasigma
            newsigmasq(1) = newsigmasq(1) * betarnd(newlambdasigma, lambdasigma - newlambdasigma);
        end

        if ~isfinite(newsigmasq(1)) || newsigmasq(1) <= 0
            invalidPathProposals = invalidPathProposals + 1;
            accept = 0;
        else
            for i = 2:T
                if ~isfinite(newsigmasq(i - 1)) || newsigmasq(i - 1) <= 0
                    validPath = false;
                    invalidPathProposals = invalidPathProposals + 1;
                    break
                end

                oldmean = rhosigma / (1 - rhosigma) * lambdasigma / musigma * sigmasq(i - 1);
                newmean = newrhosigma / (1 - newrhosigma) * newlambdasigma / newmusigma * newsigmasq(i - 1);

                if ~isfinite(newmean) || ~isfinite(oldmean) || oldmean <= 0 || newmean < 0
                    validPath = false;
                    invalidPathProposals = invalidPathProposals + 1;
                    break
                end

                if (newmean > oldmean)
                    if (newmean - oldmean) >= n_condition
                        newkappasigmasq(i - 1) = kappasigmasq(i - 1) + max(0, round(normrnd((newmean - oldmean), sqrt((newmean - oldmean)))));
                        countUsage(6) = countUsage(6) + 1;
                    else
                        newkappasigmasq(i - 1) = kappasigmasq(i - 1) + poissrnd((newmean - oldmean));
                        countUsage(5) = countUsage(5) + 1;
                    end
                else
                    p_sigmasq = newmean / oldmean;
                    np = kappasigmasq(i-1) * p_sigmasq;
                    if kappasigmasq(i-1) == 0 || p_sigmasq == 0
                        newkappasigmasq(i - 1) = 0;
                        countUsage(7) = countUsage(7) + 1;
                    elseif p_sigmasq == 1
                        newkappasigmasq(i - 1) = kappasigmasq(i-1);
                        countUsage(7) = countUsage(7) + 1;
                    elseif kappasigmasq(i-1) >= n_condition && (0.5 - p_condition) <= p_sigmasq && p_sigmasq <= (p_condition + 0.5)
                        % Normal approximation
                        normalDraw = normrnd(np, sqrt(np * (1 - p_sigmasq)));
                        countUsage(2) = countUsage(2) + 1;
                        approximateCount = round(normalDraw);
                        % Clip to the valid range [0, n]
                        newkappasigmasq(i - 1) = max(0, min(kappasigmasq(i-1), approximateCount));
                    elseif kappasigmasq(i-1) >= n_condition && p_sigmasq <= (p_condition / 2) && np <= np_threshold
                        % Poisson approximation (Poisson(np))
                        val = poissrnd(np);
                        countUsage(3) = countUsage(3) + 1;
                        val = min(val, kappasigmasq(i-1));
                        newkappasigmasq(i - 1) = max(0, val);
                    elseif kappasigmasq(i-1) >= n_condition && (1 - p_sigmasq) <= (p_condition / 2) && kappasigmasq(i-1) * (1 - p_sigmasq) <= np_threshold
                        % Poisson approximation with (1 - p)
                        % Number of failures via Poisson approximation:
                        X_poiss = poissrnd(kappasigmasq(i-1) * (1 - p_sigmasq));
                        countUsage(4) = countUsage(4) + 1;
                        % Clamp to [0, n] to respect the Binomial support:
                        X_poiss = max(0, min(kappasigmasq(i-1), X_poiss));
                        % Then the "successes" is n - X:
                        newkappasigmasq(i - 1) = kappasigmasq(i-1) - X_poiss;
                    else
                        newkappasigmasq(i - 1) = binornd(kappasigmasq(i-1), p_sigmasq);
                        countUsage(1) = countUsage(1) + 1;
                    end
                end


                oldlam = kappasigmasq(i - 1) + lambdasigma;
                oldgam = rhosigma / (1 - rhosigma) * lambdasigma / musigma + lambdasigma / musigma;
                newlam = newkappasigmasq(i - 1) + newlambdasigma;
                newgam = newrhosigma / (1 - newrhosigma) * newlambdasigma / newmusigma + newlambdasigma / newmusigma;

                newsigmasq(i) = sigmasq(i) * oldgam / newgam;

                if (newlam > oldlam)
                    newsigmasq(i) = newsigmasq(i) + gamrnd(newlam - oldlam, 1 / newgam);
                elseif newlam < oldlam
                    newsigmasq(i) = newsigmasq(i) * betarnd(newlam, oldlam - newlam);
                end

                if ~isfinite(newsigmasq(i)) || newsigmasq(i) <= 0
                    validPath = false;
                    invalidPathProposals = invalidPathProposals + 1;
                    break
                end
            end

            accept = 0;
            if validPath
                [newmeankf, newvarkf, newloglike] = kf_NGAR(datastar, targetstar, Psistar, newsigmasq, rhobeta(zstar==1));
                logaccept = newloglike - loglike + 3 * (log(newlambdasigma) - log(lambdasigma)) - 1 * (newlambdasigma - lambdasigma);
                logaccept = logaccept + log(newmusigma) - log(musigma);
                logaccept = logaccept - (1 + 0.5) * log(1 + newmusigma) + (1 + 0.5) * log(1 + musigma);
                logaccept = logaccept + log(1 / rhosigma + 1 / (1 - rhosigma)) - log(1 / newrhosigma + 1 / (1 - newrhosigma));
                logaccept = logaccept + (40*0.95 - 1)*(log(newrhosigma) - log(rhosigma)) + (40*0.05 - 1)*(log(1 - newrhosigma) - log(1 - rhosigma));

                accept = 1;
                if ( isnan(logaccept) || isinf(logaccept) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end
            end
        end
    end
    sigmasqparamaccept = sigmasqparamaccept + accept;
    sigmasqparamcount = sigmasqparamcount + 1;

    if ( rand < accept )
        lambdasigma = newlambdasigma;
        musigma = newmusigma;
        rhosigma = newrhosigma;
        sigmasq = newsigmasq;
        kappasigmasq = newkappasigmasq;
        loglike = newloglike;
        meankf = newmeankf;
        varkf = newvarkf;
    end

    if ( it < 100 )
        loglambdasigmasd = loglambdasigmasd + 1 / it^0.55 * (accept - 0.3);
        loggammasigmasqsd = loggammasigmasqsd + 1 / it^0.55 * (accept - 0.3);
        logrhosigmasd = logrhosigmasd + 1 / it^0.55 * (accept - 0.3);
    else
        logscalesigmasq = logscalesigmasq + 1 / (it - 99)^0.55 * (accept - 0.3);
    end

    x1 = log(lambdasigma);
    x2 = log(musigma);
    x3 = log(rhosigma) - log(1 - rhosigma);

    sum1sigmasq(1) = sum1sigmasq(1) + x1;
    sum1sigmasq(2) = sum1sigmasq(2) + x2;
    sum1sigmasq(3) = sum1sigmasq(3) + x3;

    sum2sigmasq(1) = sum2sigmasq(1) + x1^2;
    sum2sigmasq(2) = sum2sigmasq(2) + x1*x2;
    sum2sigmasq(3) = sum2sigmasq(3) + x2^2;
    sum2sigmasq(4) = sum2sigmasq(4) + x1*x3;
    sum2sigmasq(5) = sum2sigmasq(5) + x2*x3;
    sum2sigmasq(6) = sum2sigmasq(6) + x3^2;






    counter = 0;
    for j = 1:p
        if ( zstar(j) == 1 )
            counter = counter + 1;

            xstar = [log(lambda(j)); log(mu(j)); log(rho(j))-log(1-rho(j)); log(rhobeta(j))-log(1-rhobeta(j))];

            if ( it < start_adap )
                newxstar = xstar + [exp(loglambdasd(j)); exp(logmeansd(j)); exp(logrhosd(j)); exp(logrhobetasd(j))] .* randn(4, 1);
            else
                sxx = [sum2(1, j) sum2(2, j) sum2(4, j) sum2(7, j); sum2(2, j) sum2(3, j) sum2(5, j) sum2(8, j); sum2(4, j) sum2(5, j) sum2(6, j) sum2(9, j); sum2(7, j) sum2(8, j) sum2(9, j) sum2(10, j)];
                varstar1 = (sxx - sum1(:, j) * sum1(:, j)' / (it - start_samples)) / (it - start_samples - 1);

                [proposalChol, proposalStatus] = chol(exp(logscale(j)) * varstar1);
                if proposalStatus == 0
                    newxstar = xstar + proposalChol' * randn(4, 1);
                else
                    % If the learned covariance cannot be factored, use independent parameter steps.
                    newxstar = xstar + [exp(loglambdasd(j)); exp(logmeansd(j)); exp(logrhosd(j)); exp(logrhobetasd(j))] .* randn(4, 1);
                end
            end

            newlambda = exp(newxstar(1));
            newmu = exp(newxstar(2));
            newrho = exp(newxstar(3)) / (1 + exp(newxstar(3)));
            newrhobeta = rhobeta;
            newrhobeta(j) = exp(newxstar(4)) / (1 + exp(newxstar(4)));
            newDelta = newrho / (1 - newrho) * newlambda / newmu;

            if (~all(isfinite([newlambda, newmu, newrho, newrhobeta(j), newDelta])) || ...
                newrho <= 0 || newrhobeta(j) <= 0 || ...
                (newrhobeta(j) > limit1) || (newrho > limit1) || ...
                (newlambda     < lambda_min) || (newlambda > lambda_max ) || ...
                (newmu         < mu_min    ) || (newmu     > mu_max     ))
                accept = 0;
            else
                validPath = true;
                newPsistar = Psistar;
                newkappastar = kappastar;

                newPsistar(1, counter) = Psistar(1, counter) * (lambda(j) / newlambda) * newmu / mu(j);
                if ( newlambda > lambda(j) )
                    newPsistar(1, counter) = newPsistar(1, counter) + gamrnd(newlambda - lambda(j), newmu / newlambda);
                elseif newlambda < lambda(j)
                    newPsistar(1, counter) = newPsistar(1, counter) * betarnd(newlambda, lambda(j) - newlambda);
                end

                if ~isfinite(newPsistar(1, counter)) || newPsistar(1, counter) <= 0
                    invalidPathProposals = invalidPathProposals + 1;
                    accept = 0;
                else

                    for i = 2:T
                        if ~isfinite(newPsistar(i - 1, counter)) || newPsistar(i - 1, counter) <= 0
                            validPath = false;
                            invalidPathProposals = invalidPathProposals + 1;
                            break
                        end

                        oldmean = Delta(j) * Psistar(i - 1, counter);
                        newmean = newDelta * newPsistar(i - 1, counter);

                        if ~isfinite(newmean) || ~isfinite(oldmean) || oldmean <= 0 || newmean < 0
                            validPath = false;
                            invalidPathProposals = invalidPathProposals + 1;
                            break
                        end

                        if (newmean > oldmean)
                            if (newmean - oldmean) >= n_condition
                                newkappastar(i - 1, counter) = kappastar(i - 1, counter) + max(0, round(normrnd((newmean - oldmean), sqrt((newmean - oldmean)))));
                                countUsage(6) = countUsage(6) + 1;
                            else
                                newkappastar(i - 1, counter) = kappastar(i - 1, counter) + poissrnd((newmean - oldmean));
                                countUsage(5) = countUsage(5) + 1;
                            end
                        else
                            p_star = newmean / oldmean;
                            np = kappastar(i - 1, counter) * p_star;
                            if kappastar(i - 1, counter) == 0 || p_star == 0
                                newkappastar(i - 1, counter) = 0;
                                countUsage(7) = countUsage(7) + 1;
                            elseif p_star == 1
                                newkappastar(i - 1, counter) = kappastar(i - 1, counter);
                                countUsage(7) = countUsage(7) + 1;
                            elseif kappastar(i - 1, counter) >= n_condition && (0.5 - p_condition) <= p_star && p_star <= (p_condition + 0.5)
                                % Normal approximation
                                normalDraw = normrnd(np, sqrt(np * (1 - p_star)));
                                countUsage(2) = countUsage(2) + 1;
                                approximateCount = round(normalDraw);
                                % Clip to the valid range [0, n]
                                newkappastar(i - 1, counter) = max(0, min(kappastar(i - 1, counter), approximateCount));
                            elseif kappastar(i - 1, counter) >= n_condition && p_star <= (p_condition / 2) && np <= np_threshold
                                % Poisson approximation (Poisson(np))
                                val = poissrnd(np);
                                countUsage(3) = countUsage(3) + 1;
                                val = min(val, kappastar(i - 1, counter));
                                newkappastar(i - 1, counter) = max(0, val);
                            elseif kappastar(i - 1, counter) >= n_condition && (1 - p_star) <= (p_condition / 2) && kappastar(i - 1, counter) * (1 - p_star) <= np_threshold
                                % Poisson approximation with (1 - p)
                                % Number of failures via Poisson approximation:
                                X_star = poissrnd(kappastar(i - 1, counter) * (1 - p_star));
                                countUsage(4) = countUsage(4) + 1;
                                % Clamp to [0, n] to respect the Binomial support:
                                X_star = max(0, min(kappastar(i - 1, counter), X_star));
                                % Then the "successes" is n - X:
                                newkappastar(i - 1, counter) = kappastar(i - 1, counter) - X_star;
                            else
                                newkappastar(i - 1, counter) = binornd(kappastar(i - 1, counter), p_star);
                                countUsage(1) = countUsage(1) + 1;
                            end
                        end


                        oldlam = kappastar(i - 1, counter) + lambda(j);
                        oldgam = Delta(j) + lambda(j) / mu(j);
                        newlam = newkappastar(i - 1, counter) + newlambda;
                        newgam = newDelta + newlambda / newmu;

                        newPsistar(i, counter) = Psistar(i, counter) * oldgam / newgam;

                        if (newlam > oldlam)
                            newPsistar(i, counter) = newPsistar(i, counter) + gamrnd(newlam - oldlam, 1 / newgam);
                        elseif newlam < oldlam
                            newPsistar(i, counter) = newPsistar(i, counter) * betarnd(newlam, oldlam - newlam);
                        end

                        if ~isfinite(newPsistar(i, counter)) || newPsistar(i, counter) <= 0
                            validPath = false;
                            invalidPathProposals = invalidPathProposals + 1;
                            break
                        end
                    end

                    accept = 0;
                    if validPath
                        [newmeankf, newvarkf, newloglike] = kf_NGAR(datastar, targetstar, newPsistar, sigmasq, newrhobeta(zstar==1));
                        logaccept = newloglike - loglike;
                        if (j == 1)
                            logaccept = logaccept + 2 * log(newlambda) - 2 * log(lambda(j)) - 4 * log(0.5 + newlambda) + 4 * log(0.5 + lambda(j));
                            logaccept = logaccept + log(newmu) - log(mu(j));
                            logaccept = logaccept - (1 + 0.5) * log(1 + newmu) + (1 + 0.5) * log(1 + mu(j));
                        else
                            logaccept = logaccept + 2 * log(newlambda) - 2 * log(lambda(j)) - 4 * log(0.5 + newlambda) + 4 * log(0.5 + lambda(j));
                            logaccept = logaccept + lambdastar * (log(newmu) - log(mu(j))) - lambdastar / mustar * (newmu - mu(j));
                        end
                        logaccept = logaccept + log(1 / rho(j) + 1 / (1 - rho(j))) - log(1 / newrho + 1 / (1 - newrho));
                        logaccept = logaccept + (80*0.97 - 1)*(log(newrho) - log(rho(j))) + (80*0.03 - 1)*(log(1 - newrho) - log(1 - rho(j)));
                        logaccept = logaccept + log(1 / rhobeta(j) + 1 / (1 - rhobeta(j))) - log(1 / newrhobeta(j) + 1 / (1 - newrhobeta(j)));
                        logaccept = logaccept + (80*0.97 - 1)*(log(newrhobeta(j)) - log(rhobeta(j))) + (80*0.03 - 1)*(log(1 - newrhobeta(j)) - log(1 - rhobeta(j)));

                        accept = 1;
                        if ( isnan(logaccept) || isinf(logaccept) )
                            accept = 0;
                        elseif ( logaccept < 0 )
                            accept = exp(logaccept);
                        end
                    end
                end
            end

            Psiparamaccept(j) = Psiparamaccept(j) + accept;
            Psiparamcount(j) = Psiparamcount(j) + 1;

            if ( rand < accept )
                lambda(j) = newlambda;
                mu(j) = newmu;
                rho(j) = newrho;
                rhobeta(j) = newrhobeta(j);
                Delta(j) = newDelta;
                Psistar(:, counter) = newPsistar(:, counter);
                kappastar(:, counter) = newkappastar(:, counter);
                loglike = newloglike;
                meankf = newmeankf;
                varkf = newvarkf;

            end

            if ( it < start_adap )
                loglambdasd(j) = loglambdasd(j) + 1 / it^0.55 * (accept - 0.3);
                logmeansd(j) = logmeansd(j) + 1 / it^0.55 * (accept - 0.3);
                logrhosd(j) = logrhosd(j) + 1 / it^0.55 * (accept - 0.3);
                logrhobetasd(j) = logrhobetasd(j) + 1 / it^0.55 * (accept - 0.3);
            else
                logscale(j) = logscale(j) + 1 / (it - 99)^0.55 * (accept - 0.3);
            end
        end
    end

    if ( it >= start_samples )
        x1 = log(lambda);
        x2 = log(mu);
        x3 = log(rho) - log(1 - rho);
        x4 = log(rhobeta) - log(1 - rhobeta);

        sum1(1, :) = sum1(1, :) + x1;
        sum1(2, :) = sum1(2, :) + x2;
        sum1(3, :) = sum1(3, :) + x3;
        sum1(4, :) = sum1(4, :) + x4;

        sum2(1, :) = sum2(1, :) + x1.^2;
        sum2(2, :) = sum2(2, :) + x1 .* x2;
        sum2(3, :) = sum2(3, :) + x2.^2;
        sum2(4, :) = sum2(4, :) + x1 .* x3;
        sum2(5, :) = sum2(5, :) + x2 .* x3;
        sum2(6, :) = sum2(6, :) + x3.^2;
        sum2(7, :) = sum2(7, :) + x1 .* x4;
        sum2(8, :) = sum2(8, :) + x2 .* x4;
        sum2(9, :) = sum2(9, :) + x3 .* x4;
        sum2(10, :) = sum2(10, :) + x4.^2;
    end



    Psi(:, zstar == 1) = Psistar;
    kappa(:, zstar == 1) = kappastar;

    % Draw the selected coefficients using the backward sampler.
    newbeta = zeros(T, sum(zstar));

    checkstar = 1;
    [cholstar, check] = chol(varkf(:, :, T));
    if ( check == 0 )
        newbeta(T, :) = (meankf(:, T) + cholstar' * randn(size(cholstar, 2), 1))';
    else
        error('NGAR:StateCovariance', 'A coefficient covariance is not positive definite; no output has been returned.');
    end
    for i = (T-1):-1:1
        Gkal = diag(rhobeta(zstar == 1) .* sqrt(Psistar(i + 1, :) ./ Psistar(i, :)));
        invQ = diag(1 ./ (1 - rhobeta(zstar == 1).^2) .* 1./Psistar(i + 1, :));
        % Scale the matrix to unit diagonal, invert it, then undo the scaling.
        covScale = sqrt(diag(varkf(:, :, i)));
        invvarkf = ((varkf(:, :, i) ./ (covScale * covScale')) \ eye(numel(covScale))) ./ (covScale * covScale');
        precisionfb = invvarkf + Gkal' * invQ * Gkal;
        precisionScale = sqrt(diag(precisionfb));
        varfb = ((precisionfb ./ (precisionScale * precisionScale')) \ eye(numel(precisionScale))) ./ (precisionScale * precisionScale');
        meanfb = varfb * (invvarkf * meankf(:, i) + Gkal' * invQ * newbeta(i + 1, :)');
        [cholstar, check] = chol(varfb);
        if ( check == 0 )
            newbeta(i, :) = (meanfb + cholstar' * randn(size(cholstar, 2), 1))';
        else
            error('NGAR:StateCovariance', 'A coefficient covariance is not positive definite; no output has been returned.');
        end
    end

    if all(isfinite(newbeta(:)))
        beta(:, zstar==1) = newbeta;
    else
        error('NGAR:InvalidBeta', 'The coefficient draw is nonfinite; no output has been returned.');
    end





    % Update mustar
    newmustar = mustar * exp(mugammasd * randn);

    logaccept = (p - 1) * lambdastar * (log(mustar) - log(newmustar));
    logaccept = logaccept - lambdastar * (1 / newmustar - 1 / mustar) * sum(mu(2:end));
    logaccept = logaccept + log(newmustar) - log(mustar) - 3 * log(newmustar + mumean1b) + 3 * log(mustar + mumean1b);

    accept = 1;
    if ( isnan(logaccept) || isinf(logaccept) )
        accept = 0;
    elseif ( logaccept < 0 )
        accept = exp(logaccept);
    end


    mugammaaccept = mugammaaccept + accept;
    mugammacount = mugammacount + 1;

    if ( rand < accept )
        mustar = newmustar;
    end

    newmugammasd = mugammasd + 1 / it^0.5 * (accept - 0.3);

    if ( (newmugammasd > 10^(-3)) && (newmugammasd < 10^3) )
        mugammasd = newmugammasd;
    end




    % Update lambdastar
    newlambdastar = lambdastar * exp(vgammasd * randn);

    logaccept = (p - 1) * (newlambdastar * log(newlambdastar / mustar) - lambdastar * log(lambdastar / mustar));
    logaccept = logaccept - (p - 1) * (gammaln(newlambdastar) - gammaln(lambdastar));
    logaccept = logaccept + (newlambdastar - lambdastar) * sum(log(mu(2:end)));
    logaccept = logaccept - (newlambdastar - lambdastar) / mustar * sum(mu(2:end));
    logaccept = logaccept + log(newlambdastar) - log(lambdastar) - 1 / mulambdastar * (newlambdastar - lambdastar);

    accept = 1;
    if ( isnan(logaccept) || isinf(logaccept) )
        accept = 0;
    elseif ( logaccept < 0 )
        accept = exp(logaccept);
    end

    vgammaaccept = vgammaaccept + accept;
    vgammacount = vgammacount + 1;

    if ( rand < accept )
        lambdastar = newlambdastar;
    end

    newvgammasd = vgammasd + 1 / it^0.5 * (accept - 0.3);

    if ( (newvgammasd > 10^(-3)) && (newvgammasd < 10^3) )
        vgammasd = newvgammasd;
    end


    % Store retained draws.
    if ( (it > burnin) && (mod(it - burnin, every) == 0) )
        if saveParameters.beta
            holdbeta(:, :, (it - burnin) / every) = beta;
        end
        if saveParameters.Psi
            holdPsi(:, :, (it - burnin) / every) = Psi;
        end
        if saveParameters.sigmasq
            holdsigmasq(:, (it - burnin) / every) = sigmasq;
        end
        if saveParameters.lambda
            holdlambda(:, (it - burnin) / every) = lambda;
        end
        if saveParameters.mu
            holdmu(:, (it - burnin) / every) = mu;
        end
        if saveParameters.rho
            holdrho(:, (it - burnin) / every) = rho;
        end
        if saveParameters.rhobeta
            holdrhobeta(:, (it - burnin) / every) = rhobeta;
        end
        if saveParameters.lambdasigma
            holdlambdasigma((it - burnin) / every) = lambdasigma;
        end
        if saveParameters.musigma
            holdmusigma((it - burnin) / every) = musigma;
        end
        if saveParameters.rhosigma
            holdrhosigma((it - burnin) / every) = rhosigma;
        end
        if saveParameters.lambdastar
            holdlambdastar((it - burnin) / every) = lambdastar;
        end
        if saveParameters.mustar
            holdmustar((it - burnin) / every) = mustar;
        end
    end
    if options.PrintEvery > 0 && mod(it, options.PrintEvery) == 0
        fprintf('Iteration %d of %d\n', it, numberofiterations);
        disp(['lambdastar = ' num2str(lambdastar)]);
        disp(['mustar = ' num2str(mustar)]);
        disp(['lambda sigmasq = ' num2str(lambdasigma)]);
        disp(['mean sigmasq = ' num2str(musigma)]);
        disp(['rhosigma = ' num2str(rhosigma)]);
        disp(['lambda = ' num2str(lambda)]);
        disp(['mu = ' num2str(mu)]);
        disp(['rho beta = ' num2str(rhobeta)]);
        disp(['rho = ' num2str(rho)]);
        disp(' ')
        disp(['kappa accept = ' num2str(kappaaccept/kappacount)]);
        disp(['kappasigmasq accept = ' num2str(kappalambdasigmaccept/kappasigmasqcount)]);
        disp(['Psi param accept = (' num2str(min(Psiparamaccept./Psiparamcount)) ', ' num2str(max(Psiparamaccept./Psiparamcount)) ')']);
        disp(['sigmasq param accept = ' num2str(sigmasqparamaccept/sigmasqparamcount)]);
        disp(['mugamma accept = ' num2str(mugammaaccept/mugammacount)]);
        disp(['vgamma accept = ' num2str(vgammaaccept/vgammacount)]);
        disp(['sigmasq accept = (' num2str(min(sigmasq1accept./sigmasq1count)) ', ' num2str(max(sigmasq1accept./sigmasq1count)) ')']);
        disp(' ')
    end
end

output = struct('beta', holdbeta, 'sigmasq', holdsigmasq, 'Psi', holdPsi, 'lambda', holdlambda, 'mu', holdmu, 'rhobeta', holdrhobeta, 'rho', holdrho, 'lambdasigma', holdlambdasigma, 'musigma', holdmusigma, 'rhosigma', holdrhosigma, 'lambdastar', holdlambdastar, 'mustar', holdmustar);
output = rmfield(output, parameterNames(~cellfun(@(name) saveParameters.(name), parameterNames)));

settings = options;
settings.mumean1b = mumean1b;
settings.mulambdastar = mulambdastar;
settings.burnin = burnin;
settings.numbofits = numbofits;
settings.every = every;
settings.n_condition = n_condition;
settings.p_condition = p_condition;
settings.dataSize = [T, p];
settings.targetSize = inputTargetSize;
output.diagnostics.settings = settings;
output.diagnostics.iterations = numberofiterations;
output.diagnostics.burnin = burnin;
output.diagnostics.thinning = every;
output.diagnostics.rngInitial = initialRng;
output.diagnostics.rngFinal = rng;
output.diagnostics.countBranchNames = {'exactBinomial', 'normalBinomial', ...
    'poissonSuccesses', 'poissonFailures', 'exactPoisson', 'normalPoisson', 'deterministic'};
output.diagnostics.countBranchUsage = countUsage;
output.diagnostics.approximateDraws = sum(countUsage([2, 3, 4, 6]));
output.diagnostics.invalidPathProposals = invalidPathProposals;
output.diagnostics.elapsedSeconds = toc(runTimer);
