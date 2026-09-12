function tests = test_normal_gamma
tests = functiontests(localfunctions);
end

function testRequestedSigns(testCase)
componentFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(componentFolder, 'normal-gamma'));
testCase.assertEqual(exist(fullfile(componentFolder, 'normal-gamma', ...
    'private', 'randraw.m'), 'file'), 2, 'The bundled RANDRAW source is missing.');

rng(20260909)
n = 30;
X = randn(n, 2);
y = X * [0.8; -0.7] + 0.5 * randn(n, 1);
fit = NG(X, y, 0.5, 0.5, 1, 1, 10, 20, 1, 1, 2);

testCase.verifyGreaterThanOrEqual(fit.beta(1, :), zeros(1, 20));
testCase.verifyLessThanOrEqual(fit.beta(2, :), zeros(1, 20));
end

function testOverlappingIndicesRejected(testCase)
componentFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(componentFolder, 'normal-gamma'));
X = randn(10, 2);
y = randn(10, 1);
testCase.verifyError(@() NG(X, y, 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1), ...
    'NG:OverlappingConstraints');
end
