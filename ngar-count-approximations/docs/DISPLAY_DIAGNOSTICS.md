# Display, saved parameters and diagnostics

The final input to `NGAR` may be a scalar options struct. Existing seven- and
nine-input calls continue to work.

```matlab
options.PrintEvery = 100;
options.SaveParameters = {'beta', 'Psi'};

% Use the approximation rules.
out = NGAR(X, y, 2, 0.1, 1000, 1000, 1, 100, 0.25, options);

% Use ordinary MATLAB count samplers.
out = NGAR(X, y, 2, 0.1, 1000, 1000, 1, options);
```

## `PrintEvery`

`PrintEvery` defaults to 20 completed iterations. Set it to zero to suppress
progress output. Each display contains the iteration, current parameter values
and the acceptance summaries used in the original working code.

## `SaveParameters`

`SaveParameters` defaults to `'all'`. It accepts one parameter name, a string
vector or a cell array of character vectors. Matching is case-insensitive; the
returned fields keep the order and spelling below.

```matlab
{'beta','sigmasq','Psi','lambda','mu','rhobeta','rho', ...
 'lambdasigma','musigma','rhosigma','lambdastar','mustar'}
```

Only the requested histories are allocated and returned. All parameters are
still sampled. Use an empty cell array (`{}`) to return diagnostics only. This
option controls in-memory output and does not write a MAT file.

## Diagnostics

`out.diagnostics` is always returned.

| Field | Contents |
| --- | --- |
| `settings` | Effective options, sampler inputs, data dimensions and the original size of `target` |
| `iterations` | Total completed iterations, including burn-in and thinned-away draws |
| `burnin`, `thinning` | Effective burn-in and thinning values |
| `elapsedSeconds` | Wall-clock time including progress display |
| `rngInitial`, `rngFinal` | MATLAB random-number states before the first draw and after the run |
| `countBranchNames` | Labels corresponding to `countBranchUsage` |
| `countBranchUsage` | Number of count values produced by each count-sampling branch |
| `approximateDraws` | Total values produced by the four explicit approximation branches |
| `invalidPathProposals` | Proposed scale paths rejected by the validity checks |

The count-branch labels are `exactBinomial`, `normalBinomial`,
`poissonSuccesses`, `poissonFailures`, `exactPoisson`, `normalPoisson` and
`deterministic`. "Exact" here distinguishes ordinary MATLAB count samplers from
the explicit approximations; it is not a broader numerical-accuracy claim.
Counts include values generated inside proposals that are later rejected. They
do not include gamma, beta, geometric or coefficient-normal draws.

The input data are not duplicated in `settings`. To repeat a run, retain the
same data and software version and restore the initial state:

```matlab
rng(out.diagnostics.rngInitial)
```

The elapsed time will naturally vary.
