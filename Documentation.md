# Documentation

- [Data format](#format)
- [Reconstruction](#reconstruction)
	* [rt_dm_reconstruct](#rt_dm_reconstruct)
	* [rt_significance](#rt_significance)
- [Information and fidelity analysis](#fidelity)
	* [rt_infomatrix](#rt_infomatrix)
	* [rt_dm_theory](#rt_dm_theory)
	* [rt_dm_theory_loss](#rt_dm_theory_loss)

<a name="format"/>
## Data format

Measurement protocol is defined by a set of complementary measurement experiments over density matrix ![rho](https://latex.codecogs.com/svg.latex?%5Crho). Every experiment is repeated many times and has several possible measurement outcomes. The probability to get ![k](https://latex.codecogs.com/svg.latex?k)-th outcome is determined by the measurement operator ![M_k](https://latex.codecogs.com/svg.latex?M_k) as ![p_k=trace(rho*M_k)](https://latex.codecogs.com/svg.latex?p_k%3D%5Ctext%7BTr%7D%28%5Crho%20M_k%29). The set of measurement operators and the number of experiments repetitions define the **_measurement protocol_**. The number of observations for each outcome define the **_measurement results_**. The following code describe the required data format.
```
proto{j}(:,:,k) % Measurement operator matrix for k-th outcome in j-th experiment
nshots(j) % Number of j-th experiment repetitions
clicks{j}(k) % Number of k-th outcome observations in j-th experiment
```

<a name="reconstruction"/>
## Reconstruction

<a name="rt_dm_reconstruct"/>
### rt_dm_reconstruct

Reconstruct the quantum state density matrix by the results of `m` complementary measurements.

j-th experiment has `l_j` possible outcomes. Total number of j-th experiment runs is `nshots(1,j)`. `clicks{j}(k,1)` is the number of observed events which correspond to measurement operator `proto{j}(:,:,k)`. The following conditions must be satisfied: `length(proto) == length(clicks) == length(nshots) = m`, `size(proto{j},3) == length(click{j}) == l_j`.

- `dm = rt_dm_reconstruct(clicks,proto,nshots)` reconstructs density matrix by the set of `clicks` obtained in each measurement experiment defined by protocol `proto` and number of runs `nshots`
- `dm = rt_dm_reconstruct(clicks,proto)` reconstructs for `nshots(j) = sum(clicks{j})`
- `dm = rt_dm_reconstruct( ___ ,Name,Value)` specifies additional parameters for reconstruction. For example, `rt_dm_reconstruct(clicks,proto,'Rank',1)` reconstructs a pure (rank-1) quantum state
- `[dm,rinfo] = rt_dm_reconstruct( ___ )` also returns additional reconstruction informations

Output:
- `dm` &ndash; density matrix
- `rinfo` &ndash; structure array with the following fields:
	* `iter` &ndash; number of iteration taken to find likelihood maximum
	* `rank` &ndash; rank of the quantum state model
	* `pval`, `chi2`, `df` &ndash; p-value, chi-squared value and number of degrees of freedom of the statistical chi-squared test
	* `dm_r`, `info_r` &ndash; cell arrays of density matrices and reconstruction information for every value of rank below optimal one (define if `Rank = 'auto'`)

List of available parameters in name-value pairs:
- `Rank` (integer) &ndash; rank of the density matrix model or 'auto' for automatic rank estimation using chi-squred test: Default: `'auto'`
- `Normalize` (boolean) &ndash; normalize output density matrix to unity. Default: `true`
- `PinvOnly` (boolean) &ndash; reconstruct density matrix by the pseudo-inversion only. Default: `false`
- `SignificanceLevel` (float) &ndash; significance level for the automatic rank definition (define if `Rank = 'auto'`). Default: `0.05`
- `Alpha` (float) &ndash; regularization coefficient. Defines the weight of previous iteration result (`Cnew -> (1-Alpha)*Cnew + Alpha*Cprev`). Default: `0.5`
- `Tol` (float) &ndash; termination tolerance on state matrix. Default: `1e-8`
- `MaxIter` (integer) &ndash; maximum number of iterations. Default: `1e6`
- `Display` (bool) &ndash; display iterations. Default: `false`

<a name="rt_significance"/>
### rt_significance

Calculate density matrix statistical significance using chi-squared test

- `pval = rt_significance(dm,clicks,proto,nshots)` calculates p-value of the chi-squared test for the density matrix `dm` and measurements results specified by `clicks`, `proto` and `nshots` (see [rt_dm_reconstruct](#rt_dm_reconstruct))
- `pval = rt_significance(dm,clicks,proto)` calculates for `nshots(j) = sum(clicks{j})`
- `pval = rt_significance( ___ ,Name,Value)` specifies additional parameters for reconstruction
- `[pval,chi2,df] = rt_significance( ___ )` also returns chi-squared value and the number of degrees of freedom

Output:
- `pval` &ndash; p-value of the chi-squared test
- `chi2` &ndash; chi-squared value
- `df` &ndash; number of degrees of freedom

List of available parameters in name-value pairs:
- `FromClicks` (boolean) &ndash; indicates that `dm` was obtained from measurement data. Default: `true`
- `Rank` (integer) &ndash; `dm` model rank (define if `FromClicks = true`). Default: `rank(dm)`
- `IsProcess` (boolean) &ndash; indicates that `dm` is the ![chi](https://latex.codecogs.com/svg.latex?%5Cchi)-matrix of some quantum process (define if `FromClicks = true`). Default: `false`
- `NormalizeDM` (boolean) &ndash; normalize `dm` to make the sum of all expected counts equal to the sum of observed counts. Default: `false`
- `MeasDF` (integer) &ndash; number of measurements degrees of freedom. Default: `N-Npovm-Nnorm`, where `N` is the total number of all possible outcomes in all experiments, `Npovm` is the number of POVM-measurements, `Nnorm = 1` if `NormalizeDM = true` and `Npovm < length(proto)` and `Nnorm = 0` otherwise

<a name="fidelity"/>
## Information and fidelity analysis

<a name="rt_infomatrix"/>
### rt_infomatrix

Calculate total Fisher information matrix

- `H = rt_infomatrix(dm,proto,nshots,r)` calculates information matrix for the density matrix `dm`, measurements type specified by `proto` and `nshots` (see [rt_dm_reconstruct](#rt_dm_reconstruct)) and rank `r` quantum state model
- `H = rt_infomatrix(dm,proto,nshots)` calculates for `r = rank(dm)`

Output:
- `H` &ndash; total Fisher information matrix

<a name="rt_dm_theory"/>
### rt_dm_theory

Calculate parameters of theoretical infidelity distribution

- `d = rt_dm_theory(dm,proto,nshots,r)` calculates infidelity distribution parameters for the density matrix `dm`, measurements type specified by `proto` and `nshots` (see [rt_dm_reconstruct](#rt_dm_reconstruct)) and rank `r` quantum state model
- `d = rt_dm_theory(dm,proto,nshots)` calculates for `r = rank(dm)`

Output:
- `d` &ndash; column-vector of infidelity distribution parameters

<a name="rt_dm_theory_loss"/>
### rt_dm_theory_loss

Calculate loss function

- `loss = rt_dm_theory_loss(dm,proto,r)` calculates loss function value for the density matrix `dm`, measurements type specified by `proto` (see [rt_dm_reconstruct](#rt_dm_reconstruct)) and rank `r` quantum state model
- `loss = rt_dm_theory_loss(dm,nshots)` calculates for `r = rank(dm)`

Output:
- `loss` &ndash; loss function value