% Small synthetic example for the sign-constrained Normal-Gamma sampler.
% Requires Statistics and Machine Learning Toolbox; randraw.m is bundled.

exampleFolder = fileparts(mfilename('fullpath'));
componentFolder = fileparts(exampleFolder);
addpath(fullfile(componentFolder, 'normal-gamma'));

rng(20260909)
n = 50;
X = randn(n, 2);
y = X * [0.8; -0.7] + 0.5 * randn(n, 1);

fit = NG(X, y, 0.5, 0.5, 1, 1, 50, 100, 1, 1, 2);

assert(all(fit.beta(1, :) >= 0))
assert(all(fit.beta(2, :) <= 0))
disp('Normal-Gamma example respected the requested signs.')
