function tests = test_NGAR
%TEST_NGAR Focused regression checks for the minimal NGAR revision.
tests = functiontests(localfunctions);
end

function setupOnce(t)
t.TestData.path = path;
t.TestData.rng = rng;
componentFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(componentFolder, 'src'));
end

function teardownOnce(t)
path(t.TestData.path);
rng(t.TestData.rng);
end

function testTargetOrientationAndFields(t)
X = [ones(10,1),linspace(-1,1,10)']; y = 0.1*sin((1:10)');
rng(20); evalc('a = NGAR(X,y,2,0.1,10,6,1,100,0.25);');
rng(20); evalc('b = NGAR(X,y'',2,0.1,10,6,1,100,0.25);');
verifyEqual(t,rmfield(a,'diagnostics'),rmfield(b,'diagnostics'));
verifySize(t,a.beta,[10,2,6]);
verifySize(t,a.Psi,[10,2,6]);
verifySize(t,a.sigmasq,[10,6]);
verifyEqual(t,numel(fieldnames(a)),13);
verifyTrue(t,all(isfinite(a.beta),'all'));
verifyTrue(t,all(a.Psi>0,'all') && all(a.sigmasq>0,'all'));
verifyTrue(t,all(a.lambda>=0.1 & a.lambda<=100,'all'));
verifyTrue(t,all(a.mu>=0.01 & a.mu<=100,'all'));
end

function testDisabledApproximationsAndThinning(t)
X = [ones(8,1),linspace(-1,1,8)']; y = 0.1*cos((1:8)');
rng(21); evalc('a = NGAR(X,y,2,0.1,8,4,2);');
rng(21); evalc('b = NGAR(X,y,2,0.1,8,8,1,Inf,0.25);');
verifyEqual(t,a.beta,b.beta(:,:,2:2:8));
verifyEqual(t,a.sigmasq,b.sigmasq(:,2:2:8));
end

function testBothAdaptationStages(t)
X = [ones(8,1),linspace(-1,1,8)']; y = 0.1*sin((1:8)');
rng(22); lastwarn('');
evalc('a = NGAR(X,y,2,0.1,1050,5,1,100,0.25);');
[~, warningId] = lastwarn;
verifyEmpty(t,warningId);
verifyTrue(t,all(isfinite(a.beta),'all'));
verifySize(t,a.beta,[8,2,5]);
end

function testApproximationAndLargerBlock(t)
X = [ones(8,1),linspace(-1,1,8)']; y = 0.1*sin((1:8)');
rng(23); evalc('a = NGAR(X,y,0.1,1,120,5,1,1,0.5);');
verifyTrue(t,all(isfinite(a.beta),'all'));
rng(24); X = [ones(6,1),randn(6,9)]; y = 0.1*randn(6,1);
evalc('b = NGAR(X,y,2,0.1,15,5,1,100,0.25);');
verifySize(t,b.beta,[6,10,5]);
verifyTrue(t,all(isfinite(b.beta),'all'));
end

function testSmallestModelAndInitialBounds(t)
rng(25); a = NGAR(ones(2,1),[0;0.1],0.001,0.001,0,1,1);
verifySize(t,a.beta,[2,1]);
verifyGreaterThanOrEqual(t,a.lambda,0.1);
verifyGreaterThanOrEqual(t,a.mu,0.01);
end

function testInvalidInputsAndWarningState(t)
X = ones(3,1); y = [0;0.1;0];
verifyError(t,@() NGAR(X,y,2,0.1,0,1,1,100),'NGAR:MissingThreshold');
verifyError(t,@() NGAR(1,1,2,0.1,0,1,1),'NGAR:TooFewObservations');
verifyError(t,@() NGAR(X,y,2,0.1,0,1,0),'MATLAB:NGAR:expectedPositive');
w = warning;
rng(26); NGAR(X,y,2,0.1,0,1,1);
verifyEqual(t,warning,w);
end

function testFilterMatchesOriginalEquations(t)
rng(27);
T = 12; p = 3;
X = randn(T,p); y = randn(T,1); Psi = exp(randn(T,p));
sigma = exp(randn(T,1)); rho = [0.8,0.9,0.95];
[m,C,ll] = kf_NGAR(X,y,Psi,sigma,rho);
[m0,C0,ll0] = originalFilter(X,y,Psi,sigma,rho);
verifyEqual(t,m,m0,'AbsTol',1e-12);
verifyEqual(t,C,C0,'AbsTol',1e-12);
verifyEqual(t,ll,ll0,'AbsTol',1e-12);
% Empty selected coefficient blocks retain the observation-only likelihood.
[~,~,ll] = kf_NGAR(zeros(T,0),y,zeros(T,0),sigma,[]);
verifyEqual(t,ll,-0.5*sum(log(sigma)+y.^2./sigma),'AbsTol',1e-12);
end

function testBackwardMeanAlgebra(t)
rng(28);
for k = 1:30
    A = randn(3); C = A*A'+eye(3);
    G = diag(0.5+0.4*rand(3,1)); Qinv = diag(1+rand(3,1));
    m = randn(3,1); next = randn(3,1);
    % Baseline expression deliberately retained here as the independent check.
    V = inv(inv(C)+G'*Qinv*G);
    reference = V*(inv(C)*m+G'*Qinv*next);
    covScale = sqrt(diag(C));
    invC = ((C./(covScale*covScale'))\eye(3))./(covScale*covScale');
    precision = invC+G'*Qinv*G;
    precisionScale = sqrt(diag(precision));
    actualV = ((precision./(precisionScale*precisionScale'))\eye(3))./(precisionScale*precisionScale');
    actual = actualV*(invC*m+G'*Qinv*next);
    verifyEqual(t,actualV,V,'AbsTol',1e-12);
    verifyEqual(t,actual,reference,'AbsTol',1e-12);
end
end

function [m,C,ll] = originalFilter(X,y,Psi,sigma,rho)
% Original Kalman equations, with the scalar inverse evaluated as a reciprocal.
[T,p] = size(X); m=zeros(p,T); C=zeros(p,p,T); ll=0;
a = zeros(p,1); P=diag(Psi(1,:));
for i=1:T
    if i>1
        Q=diag((1-rho.^2).*Psi(i,:));
        G=diag(rho.*sqrt(Psi(i,:)./Psi(i-1,:)));
        a=G*m(:,i-1); P=G*C(:,:,i-1)*G'+Q;
    end
    x=X(i,:); e=y(i)-x*a; invF=1/(sigma(i)+x*P*x');
    m(:,i)=a+P*x'*invF*e;
    C(:,:,i)=P-P*x'*invF*x*P;
    ll=ll-0.5*e^2*invF+0.5*log(invF);
end
end
