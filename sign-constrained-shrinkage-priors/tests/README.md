# Software checks

After running `install_overlays.R`, run these commands from the component folder:

```sh
Rscript --vanilla tests/smoke_installed.R
Rscript --vanilla examples/example_constraints.R
```

The installer also runs `smoke_installed.R` automatically. These scripts select the local library with `use_overlays.R` and call the installed package functions. They do not replace functions by sourcing the overlay files. The test reports package versions and paths, so an ordinary installed package cannot silently stand in for a modified one.

The R checks cover:

- constrained draws from `bayesreg`, `dl()` and `r2d2marg()`;
- the compiled `Boom`/`BoomSpikeSlab` changes on the Gaussian classic SSVS route;
- overlapping and out-of-range constraint indices; and
- rejection of unsupported adaptive SSVS requests.

The test data favour slopes with the opposite signs to the requested restrictions. The spike-and-slab sign test holds all coefficients in the model so it cannot pass just because the constrained coefficients were excluded. These are short execution and interface checks. They do not assess MCMC convergence, mixing or equality with an independently implemented constrained posterior.

The installation workflow and R checks were tested with R 4.5.1 on Apple silicon macOS. All four modified packages were built from pristine pinned sources into an empty local library; dependencies already present in the system library were reused. Installation with a completely empty dependency library, and Windows/Linux builds, have not been verified.

`test_normal_gamma.m` checks the MATLAB constraints using the bundled `randraw.m`. Run it from the component root with:

```matlab
results = runtests('tests/test_normal_gamma.m');
assertSuccess(results)
```
