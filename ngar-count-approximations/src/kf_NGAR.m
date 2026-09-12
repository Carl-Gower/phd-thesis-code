function [meankf, varkf, loglike] = kf_NGAR(data, y, Psi, sigmasq, rhobeta)

%KF_NGAR Kalman filter used by the NGAR dynamic regression sampler.
%   This helper is adapted from the code accompanying Kalli and Griffin
%   (2014). Published on GitHub with written permission from Maria Kalli and
%   Jim E. Griffin. No named standard software licence is asserted; see
%   PROVENANCE.md and LICENSE.md.

p = size(data, 2);
T = length(y);

meankf=zeros(p, T);
varkf=zeros(p, p, T);

aminus = zeros(p, 1);
Pminus = diag(Psi(1, :));
datastar = data(1, :);
e = y(1) - datastar * aminus;
invF = 1 / (sigmasq(1) + datastar * Pminus * datastar');
meankf(:, 1) = aminus + Pminus * datastar' * invF * e;
varkf(:, :, 1) = Pminus - Pminus * datastar' * invF * datastar * Pminus;
loglike = - 0.5 * e^2 * invF + 0.5 * log(invF);

for i = 2:T
    datastar = data(i, :);
    Q = diag((1 - rhobeta.^2) .* Psi(i, :));
    Gkal = diag(rhobeta .* sqrt(Psi(i, :) ./ Psi(i - 1, :)));
    aminus = Gkal * meankf(:, i - 1);
    Pminus = Gkal * varkf(:, :, i - 1) * Gkal' + Q;
    e = y(i) - datastar * aminus;
    invF = 1 / (sigmasq(i) + datastar * Pminus * datastar');
    meankf(:, i) = aminus + Pminus * datastar' * invF * e;
    varkf(:, :, i) = Pminus - Pminus * datastar' * invF * datastar * Pminus;
    loglike = loglike - 0.5 * e^2 * invF + 0.5 * log(invF);
end
