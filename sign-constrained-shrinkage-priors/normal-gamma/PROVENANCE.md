# Normal-Gamma provenance

`NG.m` is derived from a MATLAB implementation supplied by Jim E. Griffin for the Normal-Gamma regression prior associated with Griffin and Brown (2010).

Carl Gower added the `pos_con` and `neg_con` interface and the coordinate-wise truncated-normal coefficient update, along with input checks and numerical safeguards. The original unconstrained block update is retained when neither argument is supplied.

The source copy supplied for this work did not include a named standard software licence. Jim E. Griffin gave Carl Gower written permission by email to put the modified code on GitHub. The correspondence is retained privately and the repository does not attempt to expand or restrict its terms. See [`LICENSE.md`](LICENSE.md).

Griffin, J. E. and Brown, P. J. (2010). Inference with normal-gamma prior distributions in regression problems. *Bayesian Analysis*, 5(1), 171–188. <https://doi.org/10.1214/10-BA507>

## RANDRAW dependency

[`private/randraw.m`](private/randraw.m) is an unmodified copy from Alex Bar-Guy's [RANDRAW submission on MATLAB Central](https://uk.mathworks.com/matlabcentral/fileexchange/7309-randraw), File Exchange release 1.4.0.0, published 6 March 2013. The source also credits Alexander Podgaetsky. The release was downloaded on 11 September 2026 from the [MathWorks source archive](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/7309/versions/5/download/zip). Its internal header uses an older version number; the release number here is the File Exchange version.

The source and its [BSD 2-Clause licence](private/LICENSE-randraw.txt) are preserved as supplied. This licence is separate from the permission for `NG.m`. Keeping the dependency in MATLAB's `private/` folder lets `NG.m` use it without a separate path setting and avoids replacing any other copy on the user's MATLAB path.
