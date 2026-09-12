function [output] = NG(data, target, lambdastar, gamma, uplambda, upgamma, burnin, numbofits, every, pos_con, neg_con)
%NG Normal-gamma Bayesian regression with optional sign constraints.
%   OUT = NG(X,Y,LAMBDASTAR,GAMMA,UPLAMBDA,UPGAMMA,BURNIN,NDRAWS,THIN)
%   retains the original unconstrained coefficient update. POS_CON and
%   NEG_CON optionally give one-based column indices of X whose coefficients
%   are restricted to be positive or negative.
%
%   Adapted by Carl Gower from MATLAB code supplied by Jim E. Griffin for the
%   normal-gamma prior of Griffin and Brown (2010). Published on GitHub with
%   written permission from Jim E. Griffin. No named standard software licence
%   is asserted; see PROVENANCE.md and LICENSE.md.
%   Requires Statistics and Machine Learning Toolbox. The BSD-licensed
%   randraw.m dependency is included in private/.

narginchk(9, 11)
if nargin < 10, pos_con = []; end
if nargin < 11, neg_con = []; end

[n,p]=size(data);
validateattributes(data, {'numeric'}, {'2d','real','finite','nonempty'}, mfilename, 'data')
validateattributes(target, {'numeric'}, {'vector','real','finite','numel',n}, mfilename, 'target')
validateattributes(lambdastar, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'lambdastar')
validateattributes(gamma, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'gamma')
validateattributes(uplambda, {'numeric','logical'}, {'scalar'}, mfilename, 'uplambda')
validateattributes(upgamma, {'numeric','logical'}, {'scalar'}, mfilename, 'upgamma')
validateattributes(burnin, {'numeric'}, {'scalar','integer','nonnegative'}, mfilename, 'burnin')
validateattributes(numbofits, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'numbofits')
validateattributes(every, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'every')
target = target(:);
pos_con = validate_constraint_indices(pos_con, p, 'pos_con');
neg_con = validate_constraint_indices(neg_con, p, 'neg_con');
if ~isempty(intersect(pos_con, neg_con))
    error('NG:OverlappingConstraints', 'pos_con and neg_con must not overlap.')
end
if exist(fullfile(fileparts(mfilename('fullpath')), 'private', 'randraw.m'), 'file') ~= 2
    error('NG:MissingRandraw', 'The bundled private/randraw.m is missing. Restore the complete normal-gamma folder.')
end

lambda = lambdastar * ones(p, 1);
a1 = 0;
b1 = 0;
lambdaaccept = 0;
lambdacount = 0;
beta = randn(p, 1);
alpha = mean(target); % Initialize alpha
sigmasq = var(target);
xtx = [ones(n, 1) data]' * [ones(n, 1) data]; % Size (p+1)x(p+1)

% Initialization for medstar
if ( p > n )
    betastar = pinv(data'*data)*data'*(target-mean(target));
    medstar = 1 / (mean(betastar.^2) * p / n);
else
    betastar = (data' * data) \ (data' * (target - mean(target)));
    medstar = 1 / mean(betastar.^2);
end

scale1 = 0;
for i = 1:p
    scale1 = scale1 + var(data(:, i));
end

Psi = lambda / gamma * 0.01;
numberofiterations = burnin + every * numbofits;

% Storage
holdPsi = zeros(p, numbofits);
holdbeta = zeros(p, numbofits);
holdsigmasq = zeros(1, numbofits);
holdgamma = zeros(1, numbofits);
holdlambda = zeros(1, numbofits);
holdalpha = zeros(1, numbofits);

gamma = lambdastar;
lambdasd = 0.01;

% -------------------------------------------------------------------------
% Main MCMC Loop
% -------------------------------------------------------------------------
for it = 1:numberofiterations
    if ( mod(it, 1000) == 0 )
        disp(['it = ' num2str(it)]);
        disp(['lambda = ' num2str(lambdastar)]);
        disp(['gamma = ' num2str(gamma)]);
        disp(['lambda accept = ' num2str(lambdaaccept/lambdacount)]);
        disp(' ');
    end
    xhat = 10^8;

    % =====================================================================
    % STEP 1: Update Beta and Alpha
    % =====================================================================

    % Retain the original block update when no signs are constrained.
    if isempty(pos_con) && isempty(neg_con)

        if ( max(Psi) / min(Psi) > xhat )
            zstar = Psi > min(Psi)*(xhat/10);
            err = target - data(:, zstar==0) * beta(zstar==0);
            selected = [1 1+find(zstar'==1)];
            precision = xtx(selected, selected) + ...
                sigmasq * diag([0; 1./Psi(zstar==1)]);
            expec = precision \ ([ones(n, 1) data(:, zstar==1)]' * err);
            varstar = precision \ eye(size(precision));
            varstar = (varstar + varstar') / 2;

            [cholstar, check] = chol(sigmasq * varstar);
            if (check == 0 )
                cholstar = cholstar';
                x = expec + cholstar * randn(sum(zstar)+1, 1);
                alpha = x(1);
                beta(zstar==1) = x(2:end);
            end
            zstar = Psi < min(Psi)*(xhat*10);
            err = target - alpha  - data(:, zstar==0) * beta(zstar==0);
            selected = 1 + find(zstar'==1);
            precision = xtx(selected, selected) + ...
                sigmasq * diag(1./Psi(zstar==1));
            expec = precision \ (data(:, zstar==1)' * err);
            varstar = precision \ eye(size(precision));
            varstar = (varstar + varstar') / 2;
            [cholstar, check] = chol(sigmasq * varstar);
            if ( check == 0 )
                cholstar = cholstar';
                x = expec + cholstar * randn(sum(zstar), 1);
                beta(zstar==1) = x;
            end
        else
            precision = xtx + sigmasq * diag([0; 1./Psi]);
            expec = precision \ ([ones(n, 1) data]' * target);
            varstar = precision \ eye(size(precision));
            varstar = (varstar + varstar') / 2;

            [cholstar, check] = chol(sigmasq * varstar);
            if ( check == 0 )
                cholstar = cholstar';
                x = expec + cholstar * randn(p + 1, 1);
                alpha = x(1);
                beta = x(2:(p+1));
            end
        end

    else
        % =================================================================
        % Coordinate-wise Gibbs update under the sign restrictions.
        % =================================================================

        % Calculate the joint posterior precision matrix.
        % Q_unscaled = X'X + diag(Prior_Precision)
        % Prior precision for alpha is 0 (flat prior), for beta_j is 1/Psi_j
        Q_prior = diag([0; 1./Psi]);
        Q_unscaled = xtx + Q_prior;

        % Solve for the unconstrained mean without forming another inverse.
        xty = [ones(n, 1) data]' * target;
        mu_vec = Q_unscaled \ xty;

        % Scaled Precision Matrix for Conditional Distributions
        Q = Q_unscaled / sigmasq;

        % Current state vector [alpha; beta]
        theta = [alpha; beta];

        % Update alpha and each regression coefficient in turn.
        for j = 1:(p+1)
            % Conditional variance and standard deviation
            % var(theta_j | theta_-j) = 1 / Q_jj
            Q_jj = Q(j,j);
            sigma_cond = sqrt(1 / Q_jj);

            % Conditional mean
            % mu_cond = mu_j - (1/Q_jj) * sum_{k!=j} Q_jk * (theta_k - mu_k)

            % Compute cross term sum_{k!=j} Q_jk * (theta_k - mu_k)
            % Efficient calculation: (Q_row_j * (theta - mu)) - Q_jj*(theta_j - mu_j)
            diff_vec = theta - mu_vec;
            cross_term = (Q(j,:) * diff_vec) - Q_jj * diff_vec(j);

            mu_cond = mu_vec(j) - (1/Q_jj) * cross_term;

            % 3. Apply Constraints
            % Map j back to beta index: index 1 is alpha, index k+1 is beta(k)
            lower = -inf;
            upper = inf;

            if j > 1
                beta_idx = j - 1;
                if ismember(beta_idx, pos_con)
                    lower = 0;
                elseif ismember(beta_idx, neg_con)
                    upper = 0;
                end
            end

            % 4. Sample
            if lower == -inf && upper == inf
                theta(j) = mu_cond + sigma_cond * randn();
            else
                theta(j) = rand_trunc_normal(mu_cond, sigma_cond, lower, upper);
            end
        end

        % Unpack theta back to alpha and beta
        alpha = theta(1);
        beta = theta(2:end);
    end

    % =====================================================================
    % STEP 2: Update Psi (Local Shrinkage)
    % =====================================================================
    badPsi = false;     % track any failed Psi updates this iteration
    for j = 1:p
        if (beta(j)^2 < 10^(-5))
            check = 0;
            if ( lambda(j) < 0.5 )
                while (check == 0)
                    Psi(j) = 1 / gamrnd(0.5-lambda(j), 1 / (0.5 * beta(j)^2));
                    check = rand < exp(- gamma * Psi(j));
                end
            else
                while (check == 0)
                    Psi(j) = gamrnd(lambda(j)-0.5, 1 / gamma);
                    check = rand < exp(- 0.5 * (beta(j)^2) / Psi(j));
                end
            end
        else
            try
                Psi(j) = randraw('gig',[lambda(j)-0.5, beta(j)^2, 2*gamma]);
            catch
                Psi(j) = 1e-8;
                badPsi = true;
            end
        end
    end
    % --- FIX: Clamp Psi to avoid infinite precision in Step 1 ---
    % 1e-15 implies a precision of 1e15, which is numerically safe
    Psi = max(Psi, 1e-15);

    mask_bad = ~isfinite(Psi); % Check for NaNs/Infs
    if any(mask_bad)
        Psi(mask_bad) = 1e-8;
        badPsi = true;
    end

    % =====================================================================
    % STEP 3: Update SigmaSq
    % =====================================================================
    sigmasq = 1 / gamrnd(a1 + 0.5 * n, 1 / (b1 + 0.5 * sum((target - alpha - data * beta).^2)));

    % =====================================================================
    % STEP 4: Update Hyperparameters (Lambda/Gamma)
    % =====================================================================
    if ( uplambda == 1 )
        if ~badPsi
            muPsi = lambdastar / gamma;
            newlambdastar = lambdastar * exp(lambdasd * randn);
            newgamma = newlambdastar / muPsi;
            newlambda = newlambdastar * ones(p, 1);
            logaccept = log(newlambdastar) - log(lambdastar) - 0.5 * (newlambdastar - lambdastar);
            logaccept = logaccept + p * newlambdastar * log(newgamma) - p * gammaln(newlambdastar);
            logaccept = logaccept - p * lambdastar * log(gamma) + p * gammaln(lambdastar);
            logaccept = logaccept + newlambdastar * sum(log(Psi)) - newgamma * sum(Psi);
            logaccept = logaccept - lambdastar * sum(log(Psi)) + gamma * sum(Psi);

            accept = 1;
            if ( logaccept < 0 )
                accept = exp(logaccept);
            end
        else
            accept = 0;
        end

            lambdasd = lambdasd + (accept - 0.3) / it;

            lambdaaccept = lambdaaccept + accept;
            lambdacount = lambdacount + 1;

            if (~badPsi) && (rand < accept)
                lambda = newlambda;
                lambdastar = newlambdastar;
                gamma = newgamma;
            end
    end
    if ( upgamma == 1 )
        gamma = gamrnd(sum(lambda) + 2, 1 / (sum(Psi) + (2 - 1) / (medstar * lambdastar)));
    end

    % Save Samples
    if ( (it > burnin) && (mod(it - burnin, every) == 0) )
        holdlambda(1, (it - burnin) / every) = lambdastar;
        holdsigmasq(:, (it - burnin) / every) = sigmasq;
        holdbeta(:, (it - burnin) / every) = beta;
        holdPsi(:, (it - burnin) / every) = Psi;
        holdgamma(1, (it - burnin) / every) = gamma;
        holdalpha(1, (it - burnin) / every) = alpha;
    end
end
output = struct('alpha', holdalpha, 'beta', holdbeta, 'sigmasq', holdsigmasq, 'Psi', holdPsi, 'lambda', holdlambda, 'gamma', holdgamma);
end

function x = rand_trunc_normal(mu, sigma, lower, upper)
    base = makedist('Normal', 'mu', mu, 'sigma', sigma);
    bounded = truncate(base, lower, upper);
    x = random(bounded);
end

function index = validate_constraint_indices(index, p, argument)
    if isempty(index)
        index = [];
        return
    end
    validateattributes(index, {'numeric'}, {'vector','real','finite','integer'}, ...
        mfilename, argument)
    index = unique(index(:)');
    if any(index < 1 | index > p)
        error('NG:InvalidConstraint', '%s indices must be between 1 and %d.', argument, p)
    end
end
