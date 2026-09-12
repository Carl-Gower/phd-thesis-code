% Small synthetic example for the NGAR sampler.

exampleFolder = fileparts(mfilename('fullpath'));
componentFolder = fileparts(exampleFolder);
addpath(fullfile(componentFolder, 'src'));

rng(42)

T = 60;
x = randn(T, 2);
X = [ones(T, 1), x];

betaTrue = [0.5 * ones(T, 1), ...
    linspace(-0.6, 0.6, T)', ...
    0.35 * sin(linspace(0, 2*pi, T)')];
y = sum(X .* betaTrue, 2) + 0.25 * randn(T, 1);

options.PrintEvery = 0;
options.SaveParameters = {'beta'};

% These short run lengths demonstrate the interface, not convergence.
fit = NGAR(X, y, 2, 0.1, 250, 100, 1, options);
betaMean = mean(fit.beta, 3);

figure
for j = 1:size(X, 2)
    subplot(size(X, 2), 1, j)
    plot(betaTrue(:, j), '--', 'DisplayName', 'Simulated')
    hold on
    plot(betaMean(:, j), 'DisplayName', 'Posterior mean')
    ylabel(sprintf('beta_%d', j - 1))
    if j == 1
        legend('Location', 'best')
    end
end
xlabel('Time')
