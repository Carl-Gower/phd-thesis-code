function tests = test_NGAR_options
%TEST_NGAR_OPTIONS Check that reporting and storage do not alter sampling.
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

function testSelectedOutputsAndRng(t)
X=[ones(8,1),linspace(-1,1,8)']; y=0.1*sin((1:8)');
opts.PrintEvery=0;
rng(51); initial=rng;
text=evalc('a=NGAR(X,y,2,0.1,10,5,1,100,0.25,opts);');
verifyEmpty(t,text);
opts.SaveParameters={'beta','Psi'};
rng(51); b=NGAR(X,y,2,0.1,10,5,1,100,0.25,opts);
verifyEqual(t,fieldnames(b),{'beta';'Psi';'diagnostics'});
verifyEqual(t,b.beta,a.beta);
verifyEqual(t,b.Psi,a.Psi);
verifyEqual(t,b.diagnostics.rngInitial,initial);
verifyEqual(t,b.diagnostics.rngFinal,a.diagnostics.rngFinal);
verifyEqual(t,b.diagnostics.countBranchUsage,a.diagnostics.countBranchUsage);
verifyEqual(t,b.diagnostics.settings.SaveParameters,{'beta','Psi'});
verifyEqual(t,b.diagnostics.settings.dataSize,[8,2]);
verifyEqual(t,b.diagnostics.settings.targetSize,[8,1]);
verifyEqual(t,b.diagnostics.settings.mumean1b,2);
verifyEqual(t,b.diagnostics.settings.numbofits,5);
verifyEqual(t,b.diagnostics.iterations,15);
verifyGreaterThan(t,b.diagnostics.elapsedSeconds,0);
verifyFalse(t,isfield(b.diagnostics.settings,'data'));
% The initial RNG state allows the same run to be reproduced.
rng(b.diagnostics.rngInitial); c=NGAR(X,y,2,0.1,10,5,1,100,0.25,opts);
verifyEqual(t,c.beta,b.beta);
end

function testOrdinaryOptionsShorthandAndDiagnosticsOnly(t)
X=ones(4,1); y=[0;0.1;0;-0.1];opts.PrintEvery=0;
opts.SaveParameters={};
rng(52);a=NGAR(X,y,2,0.1,2,3,2,opts);
rng(52);b=NGAR(X,y,2,0.1,2,3,2,Inf,0.25,opts);
verifyEqual(t,fieldnames(a),{'diagnostics'});
verifyEqual(t,a.diagnostics.rngFinal,b.diagnostics.rngFinal);
verifyEqual(t,a.diagnostics.countBranchUsage,b.diagnostics.countBranchUsage);
verifyEqual(t,a.diagnostics.approximateDraws,0);
verifyEqual(t,a.diagnostics.iterations,8);
verifyEqual(t,a.diagnostics.burnin,2);
verifyEqual(t,a.diagnostics.thinning,2);
verifyGreaterThanOrEqual(t,a.diagnostics.countBranchUsage(5),4);
end

function testDisplayIntervalDoesNotChangeDraws(t)
X=ones(4,1); y=[0;0.1;0;-0.1];opts.PrintEvery=2;
rng(53); text=evalc('a=NGAR(X,y,2,0.1,0,5,1,opts);');
verifyEqual(t,regexp(text,'Iteration \d+ of 5','match'),{'Iteration 2 of 5','Iteration 4 of 5'});
verifyFalse(t,contains(text,'NaN'));
opts.PrintEvery=0;
rng(53);b=NGAR(X,y,2,0.1,0,5,1,opts);
verifyEqual(t,rmfield(a,'diagnostics'),rmfield(b,'diagnostics'));
verifyEqual(t,a.diagnostics.rngFinal,b.diagnostics.rngFinal);
% Printing after the first completed iteration also has valid counters.
opts.PrintEvery=1;
rng(53);text=evalc('NGAR(X,y,2,0.1,0,1,1,opts);');
verifyTrue(t,contains(text,'Iteration 1 of 1'));
verifyFalse(t,contains(text,'NaN'));
end

function testParameterNamesAndValidation(t)
X=ones(3,1); y=[0;0.1;0];opts.PrintEvery=0;
opts.SaveParameters=["psi","BETA","beta"];
rng(54);a=NGAR(X,y,2,0.1,0,1,1,opts);
verifyEqual(t,fieldnames(a),{'beta';'Psi';'diagnostics'});
opts.SaveParameters='all';
rng(54);a=NGAR(X,y,2,0.1,0,1,1,opts);
verifyEqual(t,numel(fieldnames(a)),13);
state=rng;
opts.SaveParameters={'unknown'};
verifyError(t,@()NGAR(X,y,2,0.1,0,1,1,opts),'NGAR:InvalidSaveParameters');
verifyEqual(t,rng,state);
verifyError(t,@()NGAR(X,y,2,0.1,0,1,1,struct('Unknown',1)),'NGAR:UnknownOption');
verifyError(t,@()NGAR(X,y,2,0.1,0,1,1,struct('PrintEvery',-1)),'MATLAB:NGAR:expectedNonnegative');
end

function testApproximationCounters(t)
X=[ones(8,1),linspace(-1,1,8)']; y=0.1*sin((1:8)');
opts.PrintEvery=0;opts.SaveParameters={};
rng(31);a=NGAR(X,y,0.1,1,250,5,1,1,0.25,opts);
u=a.diagnostics.countBranchUsage;
verifySize(t,u,[1,7]);
verifyTrue(t,all(u>=0 & u==fix(u)));
verifyTrue(t,all(u([2,3,4,6])>0));
verifyEqual(t,a.diagnostics.approximateDraws,sum(u([2,3,4,6])));
verifyGreaterThanOrEqual(t,u(5),numel(X));
verifyLessThanOrEqual(t,sum(u),numel(X)+255*3*7);
verifyEqual(t,a.diagnostics.countBranchNames,{'exactBinomial','normalBinomial', ...
    'poissonSuccesses','poissonFailures','exactPoisson','normalPoisson','deterministic'});
end
