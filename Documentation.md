# Documentation

- [Definitions](#definitions)
- [Algorithms](#algorithms)
- [Data format](#format)
- Basic functions
	* [rt_startup](#rt_startup)
	* [rt_experiment](#rt_experiment)
	* [rt_dm_reconstruct](#rt_dm_reconstruct)
	* [rt_chi_reconstruct](#rt_chi_reconstruct)
- Fisher information and fidelity bound
	* [rt_fidelity](#rt_fidelity)
	* [rt_infomatrix](#rt_infomatrix)
	* [rt_bound](#rt_bound)
	* [rt_lossfun](#rt_lossfun)
	* [rt_gchi2pdf](#rt_gchi2pdf)
	* [rt_gchi2cdf](#rt_gchi2cdf)
	* [rt_gchi2inv](#rt_gchi2inv)
- Random entities generator
	* [rt_randstate](#rt_randstate)
	* [rt_randunitary](#rt_randunitary)
	* [rt_randprocess](#rt_randprocess)
- Protocols generators
	* [rt_proto_measurement](#rt_proto_measurement)
	* [rt_proto_preparation](#rt_proto_preparation)
	* [rt_proto_process](#rt_proto_process)
	* [rt_nshots_divide](#rt_nshots_divide)
- Utilities
	* [rt_meas_matrix](#rt_meas_matrix)
	* [rt_kron3d](#rt_kron3d)
	* [rt_pinv](#rt_pinv)
	* [rt_purify](#rt_purify)
	* [rt_prttrace](#rt_prttrace)
	* [rt_to_simplex](#rt_to_simplex)
	* [rt_pval](#rt_pval)
	* [rt_optimizer](#rt_optimizer)
- [References](#references)

## <a name="definitions">Definitions</a>
Consider a quantum state in the Hilbert space of dimension <img alt="d" src="https://render.githubusercontent.com/render/math?math=d" />. The root approach to quantum state tomography implies reconstructing a purified quantum state <img alt="psi" src="https://render.githubusercontent.com/render/math?math=\psi" /> of size <img alt="d-by-r" src="https://render.githubusercontent.com/render/math?math=d\times r" /> instead of corresponding rank-<img alt="r" src="https://render.githubusercontent.com/render/math?math=r" /> **density matrix** <img alt="rho=psi*psi" src="https://render.githubusercontent.com/render/math?math=\rho=\psi\psi^\dagger" />. Thus, matrix <img alt="psi" src="https://render.githubusercontent.com/render/math?math=\psi" /> defines a **square root** of the density matrix.

A quantum process is described by a Choi-Jamiolkowski **process matrix** <img alt="chi" src="https://render.githubusercontent.com/render/math?math=\chi" /> that is the <img alt="d^2-by-d^2" src="https://render.githubusercontent.com/render/math?math=d^2\times d^2" /> density matrix normalized by <img alt="d" src="https://render.githubusercontent.com/render/math?math=d" />. Its square root <img alt="e" src="https://render.githubusercontent.com/render/math?math=e" /> is used for the quantum process reconstruction. We also impose the **trace preserving** condition:
<p align="center"><img alt="TP condition" src="https://render.githubusercontent.com/render/math?math=\displaystyle \sum_{k=1}^{r}{E_k^\dagger E_k} = I," /></p>

where <img alt="E_k" src="https://render.githubusercontent.com/render/math?math=E_k" /> are the process **Kraus operators** and <img alt="I" src="https://render.githubusercontent.com/render/math?math=I" /> is the identity matrix.

We measure the reconstruction accuracy by Uhlmann's **fidelity** between the true state <img alt="rho" src="https://render.githubusercontent.com/render/math?math=\rho" /> and the reconstructed state <img alt="sigma" src="https://render.githubusercontent.com/render/math?math=\sigma" />:
<p align="center"><img alt="fidelity" src="https://render.githubusercontent.com/render/math?math=\displaystyle F(\rho,\sigma)=\left[\textrm{Tr}\sqrt{\sqrt\rho\sigma\sqrt\rho}\right]^2." /></p>

For a quantum process we estimate the fidelity between the process matrices normalized by unity.

According to the quantum state estimation theory the **infidelity distribution** is bounded by the general chi-squared distribution with <img alt="nu" src="https://render.githubusercontent.com/render/math?math=\nu" /> degrees of freedom [[1]](#ref1):
<p align="center"><img alt="infidelity as random variable" src="https://render.githubusercontent.com/render/math?math=\displaystyle 1-F\sim\sum_{j=1}^{\nu}{d_j\xi_j^2}," /></p>

where <img alt="d_j" src="https://render.githubusercontent.com/render/math?math=d_j" /> are positive parameters and <img alt="xi_j" src="https://render.githubusercontent.com/render/math?math=\xi_j" /> are independent random variables with standard normal distribution. The expected value and variance of infidelity are thus
<p align="center"><img alt="infidelity mean and variance" src="https://render.githubusercontent.com/render/math?math=\displaystyle \langle1-F\rangle=\sum_{j=1}^{\nu}{d_j}, \qquad \sigma_{1-F}^2=2\sum_{j=1}^{\nu}{d_j^2}." /></p>

As the infidelity lower bound is inverse proportional to the total sample size <img alt="N" src="https://render.githubusercontent.com/render/math?math=N" /> over all measurements, we also use the so-called **loss function** <img alt="L=N*mean(1-F)" src="https://render.githubusercontent.com/render/math?math=L=N\langle1-F\rangle" /> independent of the sample size.

## <a name="algorithms">Algorithms</a>
We use the _maximum likelihood_ parameters estimation (MLE). In the case of quantum process tomography MLE results in the _likelihood equation_ [[1]](#ref1):
<p align="center"><img alt="likelihood equation" src="https://render.githubusercontent.com/render/math?math=\displaystyle I\psi=J(\psi)\psi," /></p>

where
<p align="center"><img alt="I and J definitions" src="https://render.githubusercontent.com/render/math?math=\displaystyle I=\sum_j{n_jP_j},\qquad J(\psi)=\sum_j{\frac{k_j}{\textrm{Tr}(\psi\psi^\dagger P_j)}P_j}." /></p>

The sums are taken over all measurements operators in all measurements schemes. A solution is obtained by the fixed-point iteration method:
<p align="center"><img alt="fixed-point iterations" src="https://render.githubusercontent.com/render/math?math=\displaystyle \psi_{i%2B1}=(1-\mu)I^{-1}J(\psi_i)\psi_i %2B \mu\psi_i." /></p>

Here <img alt="mu" src="https://render.githubusercontent.com/render/math?math=\mu" /> is the regularization parameter. We use Moore-Penrose pseudo-inversion to get <img alt="psi_0" src="https://render.githubusercontent.com/render/math?math=\psi_0" />.

For the quantum process tomography we use the _non-convex proximal gradient ascend_ method to ensure the trace preserving constraint (for a non trace preserving processes we use the fixed-point iteration method). In particular, we consider the **Algorithm 1** from Ref. [[2]](#ref2) with Lipschitz constant equal to <img alt="N" src="https://render.githubusercontent.com/render/math?math=N" /> by default. To implement the proximal operator we make use of the fact that for a trace preserving process the block matrix of Kraus operators
<p align="center"><img alt="fixed-point iterations" src="https://render.githubusercontent.com/render/math?math=\displaystyle Q=\begin{pmatrix}E_1\\E_2\\\vdots\end{pmatrix}." /></p>

should be orthogonal. Note that this matrix is just the re-ordered version of the process matrix square root <img alt="e" src="https://render.githubusercontent.com/render/math?math=e" />. Thus, the projection of <img alt="e" src="https://render.githubusercontent.com/render/math?math=e" /> is done by the Frobenius norm projection of the corresponding non orthogonal <img alt="Q" src="https://render.githubusercontent.com/render/math?math=Q" /> onto the set of orthogonal matrices. To do this we equate all its singular values to 1.

The root approach quantum tomography implies setting the model rank <img alt="r" src="https://render.githubusercontent.com/render/math?math=r" />. If the rank is unknown we estimate it using the _adequacy criterion_ [[1]](#ref1). To do this we vary <img alt="r" src="https://render.githubusercontent.com/render/math?math=r" /> from 1 to its maximal value until the reconstruction result becomes statistically significant at some pre-chosen significance level. The procedure is also terminated if the p-value of the rank-<img alt="(r+1)" src="https://render.githubusercontent.com/render/math?math=(r%2B1)" /> model is lower than p-value of the rank-<img alt="r" src="https://render.githubusercontent.com/render/math?math=r" /> model.

## <a name="format">Data format</a>
For the **quantum state tomography** one must specify a set of complementary measurement experiments over a quantum state density matrix <img alt="rho" src="https://render.githubusercontent.com/render/math?math=\rho" />. Every experiment may be repeated many times with some sets of possible measurement outcomes. The probability to get _k_-th outcome is determined by the measurement operator <img alt="P_k" src="https://render.githubusercontent.com/render/math?math=P_k" /> as <img alt="p_k=trace(rho*P_k)" src="https://render.githubusercontent.com/render/math?math=p_k=\textrm{Tr}(\rho P_k)" />. The set of measurement operators and the number of experiments repetitions define the **_measurements protocol_**. The number of observations for each outcome define the **_measurements results_**. The following code describe the required data format.
```
proto{j}(:,:,k) % Measurement operator matrix for k-th outcome in j-th measurement scheme
nshots(j) % Number of j-th scheme repetitions
clicks{j}(k) % Number of k-th outcome observations in j-th scheme
```
One can specify `nshots` as a single integer describing the total sample size. Then it is automatically divided equally over all measurement schemes.

One can also pass a 3D array for measurements operators and a 1D array for measurements results implying only a single possible outcome in each measurement:
```
proto(:,:,j) % Measurement operator matrix for j-th measurement scheme
nshots(j) % Number of j-th scheme repetitions
clicks(j) % Number of observations in j-th scheme
```

The measurements protocol for the **quantum process tomography** is defined in a very similar way. One option is to define a 2-by-M cell array, where the first line defines the density matrix in the quantum process input. The second line defines the measurement operators in corresponding measurement experiment.
```
proto{1,j} % Input density matrix in j-th measurement scheme
proto{2,j}(:,:,k) % Output measurement operator matrix for k-th outcome in j-th scheme
nshots(j) % Number of j-th scheme repetitions
clicks{j}(k) % Number of k-th outcome observations in j-th scheme
```

Other option is to define measurement operators for the process matrix <img alt="chi" src="https://render.githubusercontent.com/render/math?math=\chi" /> and use the same data format as for quantum state tomography.

## <a name="rt_startup">rt_startup</a>
Includes the library directories in the search paths. Run before using the library.

## <a name="rt_experiment">rt_experiment</a>
The class for working with the quantum tomography data.

### Creating an object
- `ex = rt_experiment(dim, 'state')` creates an object for a quantum state tomography experiment for the Hilbert space of dimension `dim`
- `ex = rt_experiment(dim, 'process')` creates an object for a quantum process tomography experiment for the Hilbert space of dimension `dim`
- `ex = rt_experiment( ___ , stat_type)` specifies the measurements statistics type `stat_type` (default: `'auto'`). Possible values: `'poly'` (polynomial/multinomial statistics), `'poiss'` (Poisson statistics), `'auto'` (tries to automatically determine the statistics type after specifying the measurements protocol).

### Fields
- `dim` &ndash; Hilbert space dimension
- `obj_type` &ndash; object type (`'state'` or `'process'`)
- `stat_type` &ndash; measurements statistics type
- `proto` &ndash; measurements operators
- `nshots` &ndash; measurements repetitions
- `clicks` &ndash; number of observed measurements outcomes
- `vec_proto` &ndash; matrix form of the whole measurements protocol
- `vec_nshots` &ndash; matrix form of the measurements repetitions
- `vec_clicks` &ndash; matrix form of the number of observed measurements outcomes

### Methods
#### set_data
Sets the experiment data in the form of Name-Value pairs, where Name specifies the object field codename and the Value specifies the field value. The method also performs the data processing and verifying (see [Data format](#format) section).

- `ex.set_data(code1, value1, code2, value2, ___ )`

#### get_field
Returns the field data by its codename. When extracting the `vec_` field for the first time, the method performs its calculation and storing for future usage.

- `data = ex.get_field(code)` returns the field data by its codename `code`

#### get_probs_dm
Returns the vector of all-experiments outcome probabilities for the input density matrix. The method also works for a quantum process described by means of the process matrix (see [Definitions](#definitions) section).

- `p = ex.get_probs_dm(dm, tol)` returns the probabilities vector `p` by the density matrix `dm` with the tolerance `tol` (default: `0`). Values of `p` lower than `tol` are set to `tol`.

#### get_probs_sq
Returns the vector of all-experiments outcome probabilities for the input square root of the density matrix. The method also works for a quantum process described by means of the square root of the process matrix (see [Definitions](#definitions) section).

- `p = ex.get_probs_sq(sq, tol)` returns the probabilities vector `p` by the square root `sq` of a density matrix with the tolerance `tol` (default: `0`). Values of `p` lower than `tol` are set to `tol`.

#### simulate
Simulates the tomography experiment for the input density matrix. The method also works for a quantum process described by means of the process matrix (see [Definitions](#definitions) section). The simulation is performed by the pseudo-random number generation according to the `stat_type` field.

- `clicks = ex.simulate(dm)` returns the measurements outcomes cell array `clicks` (see [Data format](#format) section) by the density matrix `dm`

#### sample
Generate a pseudo-random sample according to the sample size, input probabilities and the `stat_type` field.

- `k = ex.sample(p, n)` returns the vector `k` of observed events by the vector `p` of corresponding probabilities and the sample size `n`

#### get_logL_dm
Calculates the logarithmic likelihood function for the input density matrix. The method also works for a quantum process described by means of the process matrix (see [Definitions](#definitions) section). The calculation is based on the data fields `proto`, `nshots` and `clicks` as well as on the statistics type field `stat_type`.

- `f = ex.get_logL_dm(dm)` returns the value `f` of the logarithmic likelihood function by the density matrix `dm`

#### get_logL_sq
Calculates the logarithmic likelihood function for the input square root of the density matrix. The method also works for a quantum process described by means of the square root of the process matrix (see [Definitions](#definitions) section). The calculation is based on the data fields `proto`, `nshots` and `clicks` as well as on the statistics type field `stat_type`.

- `f = ex.get_logL_sq(sq)` returns the value `f` of the logarithmic likelihood function by the square root `sq` of a density matrix

#### get_dlogL_sq
Calculates the gradient of the logarithmic likelihood function for the input square root of the density matrix. The method also works for a quantum process described by means of the square root of the process matrix (see [Definitions](#definitions) section). The calculation is based on the data fields `proto`, `nshots` and `clicks` as well as on the statistics type field `stat_type`.

- `df = ex.get_dlogL_sq(sq)` returns the matrix `df` of the logarithmic likelihood function gradient by the square root `sq` of a density matrix

#### get_chi2_dm
Calculates the chi-squared value (see [chi-squared test](https://en.wikipedia.org/wiki/Chi-squared_test)) for the input density matrix. The method also works for a quantum process described by means of the process matrix (see [Definitions](#definitions) section). The calculation is based on the data fields `proto`, `nshots` and `clicks`.

- `f = ex.get_chi2_dm(dm)` returns the value `f` of the chi-squared by the density matrix `dm`

#### get_df
Calculates the number of degrees of freedom in the [chi-squared test](https://en.wikipedia.org/wiki/Chi-squared_test) taking the rank of the model into account. The value differs for quantum state and quantum process reconstruction (specified by the `obj_type` field) and for different types of measurements statistics (`stat_type` field).

- `df = ex.get_df(obj_rank)` returns the number `df` of degrees of freedom for the model rank `obj_rank`

## <a name="rt_dm_reconstruct">rt_dm_reconstruct</a>
Reconstruct the quantum state density matrix by the results of a set of complementary measurements. The reconstruction is based on solving the likelihood equation for the density matrix square root (see [Definitions](#definitions) section) using the fixed-point iteration method (see [Algorithms](#algorithms) section).

By default the function uses the procedure of automatic rank estimation using the adequacy criterion. It probes different values of the density matrix rank from 1 to the Hilbert space dimension until the model becomes statistically significant at some significance level. The procedure is also terminated when the maximum p-value is obtained (see [Algorithms](#algorithms) section).

#### Usage
- `dm = rt_dm_reconstruct(dim, clicks, proto, nshots)` reconstructs the quantum state in the Hilbert space of dimension `dim` in the form of the density matrix `dm`. The adequacy criterion at the significance level 5% (See [Pearson's chi-squared test](https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test)) is used to estimate the state rank. Parameters `clicks`, `proto` and `nshots` describe the measurements (see [Data format](#format) section)
- `dm = rt_dm_reconstruct(dim, clicks, proto)` performs reconstruction for `nshots(j) = sum(clicks{j})`
- `dm = rt_dm_reconstruct( ___ , 'Rank', r)` reconstructs the state using the rank-`r` model
- `dm = rt_dm_reconstruct( ___ , 'Rank', 'full')` reconstructs the state using the full-rank model (`r = dim`)
- `dm = rt_dm_reconstruct( ___ , Name, Value)` specifies additional parameters for reconstruction
- `[dm, rinfo] = rt_dm_reconstruct( ___ )` also returns additional reconstruction information

#### Name-Value parameters
- `StatType` &ndash; type of the measurements statistics: `'poly'` for polynomial/multinomial statistics, `'poiss'` for Poisson statistics, and `'auto'` (default) to try to determine statistics type from the measurements operators
- `Rank` &ndash; rank of the quantum state model, or `'full'` for a full-rank model, or `'auto'` (default) for the automatic rank estimation using adequacy criterion
- `SignificanceLevel` &ndash; significance level for the automatic rank estimation (default: `0.05`)
- `GetStats` &ndash; return the chi-squared test information (default: `true` if `Rank = 'auto'`; `false` otherwise)
- `Init` &ndash; zero approximation density matrix, or `'pinv'` (default) to use the Moore-Penrose pseudo-inverse to get the zero approximation
- `RegCoeff` &ndash; regularization coefficient (default: `0.5`)
- `Tol` &ndash; optimization tolerance
- `MaxIter` &ndash; maximal number of iterations (default `1e6`)
- `Display` &ndash; display the reconstruction information (default: `false`); setting the integer value controls the frequency of the iteration procedure display update

#### Output
- `dm` &ndash; density matrix
- `rinfo` &ndash; structure array with the reconstruction information
	* `rinfo.iter` &ndash; number of optimization iteration
	* `rinfo.rank` &ndash; the reconstruction model rank
	* `rinfo.experiment` &ndash; the tomography experiment object handle (see [rt_experiment](#rt_experiment) class)
	* `rinfo.optimizer` &ndash; the optimizer object handle (see [rt_optimizer](#rt_optimizer) class)
	* `rinfo.chi2`, `rinfo.df`, `rinfo.pval` &ndash; chi-squared value, the number of model degrees of freedom, and p-value of the statistical chi-squared test (defined if `GetStats = true`)
	* `rinfo.data_r` &ndash; a cell array with the results obtained by probing models of different ranks until the model is statistically significant (defined if `Rank = 'auto'`)

## <a name="rt_chi_reconstruct">rt_chi_reconstruct</a>
Reconstruct the quantum process matrix (see [Definitions](#definitions) section) by the results of a set of complementary measurements. For the trace preserving processes the reconstruction is performed by finding the local maximum of log-likelihood function for the process matrix square root using non-convex accelerated proximal gradient ascend. For the non trace preserving processes the reconstruction is based on solving the likelihood equation using the fixed-point iteration method (see [Algorithms](#algorithms) section).

By default the function uses the procedure of automatic rank estimation using the adequacy criterion. It probes different values of the process matrix rank from 1 to the squared Hilbert space dimension until the model becomes statistically significant at some significance level. The procedure is also terminated when the maximum p-value is obtained (see [Algorithms](#algorithms) section).

#### Usage
- `chi = rt_chi_reconstruct(dim, clicks, proto, nshots)` reconstructs the quantum process in the Hilbert space of dimension `dim` in the form of the process matrix `chi`. The adequacy criterion at the significance level 5% (See [Pearson's chi-squared test](https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test)) is used to estimate the state rank. Parameters `clicks`, `proto` and `nshots` describe the measurements (see [Data format](#format) section)
- `chi = rt_chi_reconstruct(dim, clicks, proto)` performs reconstruction for `nshots(j) = sum(clicks{j})`
- `chi = rt_chi_reconstruct( ___ , 'Rank', r)` reconstructs the process using the rank-`r` model
- `chi = rt_chi_reconstruct( ___ , 'Rank', 'full')` reconstructs the process using the full-rank model (`r = dim^2`)
- `chi = rt_chi_reconstruct( ___ , Name, Value)` specifies additional parameters for reconstruction
- `[chi, rinfo] = rt_chi_reconstruct( ___ )` also returns additional reconstruction information

#### Name-Value parameters
- `TracePreserving` &ndash; is trace preserving process (default: `true`)
- `StatType` &ndash; type of the measurements statistics: `'poly'` for polynomial/multinomial statistics, `'poiss'` for Poisson statistics, and `'auto'` (default) to try to determine statistics type from the measurements operators
- `Rank` &ndash; rank of the quantum process model, or `'full'` for a full-rank model, or `'auto'` (default) for the automatic rank estimation using adequacy criterion
- `SignificanceLevel` &ndash; significance level for the automatic rank estimation (default: `0.05`)
- `GetStats` &ndash; return the chi-squared test information (default: `true` if `Rank = 'auto'`; `false` otherwise)
- `Init` &ndash; zero approximation process matrix, or `'pinv'` (default) to use the Moore-Penrose pseudoinverse to get the zero approximation
- `LipschitzConstant` &ndash; Lipschitz constant for optimization procedure, or `'ntot'` (default) to take `sum(nshots)` that usually works well; increase in case of diverging
- `Tol` &ndash; optimization tolerance
- `MaxIter` &ndash; maximal number of iterations (default `1e6`)
- `Display` &ndash; display the reconstruction information (default: `false`); setting the integer value controls the frequency of the iteration procedure display update

#### Output
- `chi` &ndash; process matrix
- `rinfo` &ndash; structure array with the reconstruction information
	* `rinfo.iter` &ndash; number of optimization iteration
	* `rinfo.rank` &ndash; the reconstruction model rank
	* `rinfo.experiment` &ndash; the tomography experiment object handle (see [rt_experiment](#rt_experiment) class)
	* `rinfo.optimizer` &ndash; the optimizer object handle (see [rt_optimizer](#rt_optimizer) class)
	* `rinfo.chi2`, `rinfo.df`, `rinfo.pval` &ndash; chi-squared value, the number of model degrees of freedom, and p-value of the statistical chi-squared test (defined if `GetStats = true`)
	* `rinfo.data_r` &ndash; a cell array with the results obtained by probing models of different ranks until the model is statistically significant (defined if `Rank = 'auto'`)

## <a name="rt_fidelity">rt_fidelity</a>
Calculates the Uhlmann's fidelity between quantum states (see [Definitions](#definitions) section).

#### Usage
- `f = rt_fidelity(dm1, dm2)` returns the fidelity between quantum states represented by density matrices `dm1` and `dm2`; the matrices are normalized by unity before computation

#### Output
- `f` &ndash; fidelity

## <a name="rt_infomatrix">rt_infomatrix</a>
Calculates the complete Fisher information matrix by the density matrix and the measurements protocol. The parameters space is the Euclidean space of real parameters of the purified quantum state. Thus, the information matrix dimension is <img alt="2dr-by-2dr" src="https://render.githubusercontent.com/render/math?math=2dr\times2dr" />. It has at least <img alt="r^2" src="https://render.githubusercontent.com/render/math?math=r^2" /> zero eigenvalues corresponded to the phase uncertainty that cannot be detected by the state measurements.

In the case of a quantum process the information matrix has dimension <img alt="2d^2r-by-2d^2r" src="https://render.githubusercontent.com/render/math?math=2d^2r\times2d^2r" /> and contains at least <img alt="r^2" src="https://render.githubusercontent.com/render/math?math=r^2" /> zero eigenvalues.

#### Usage
- `H = rt_infomatrix(dm, proto, nshots, 'state')` returns the complete Fisher information matrix `H` for the quantum state density matrix `dm` and measurements protocol described by `proto` and `nshots` (see [Data format](#format) section); the purification rank is `rank(dm)`
- `H = rt_infomatrix(chi, proto, nshots, 'process')` returns the complete Fisher information matrix for the quantum process matrix `chi`
- `H = rt_infomatrix( ___ , 'Rank', r)` specifies the purification rank `r`

#### Output
- `H` &ndash; complete Fisher information matrix

### <a name="rt_bound">rt_bound</a>
Calculates the lower bound for the variances <img alt="d_j" src="https://render.githubusercontent.com/render/math?math=d_j" /> of quantum state or quantum process parameters estimator. The parameters space is the subspace of the Euclidean space of real parameters of the purified quantum state. The subspace is orthogonal to the vectors imposed by the phase insensitivity and normalization conditions, as well as by the trace preserving condition for a quantum process.

The calculated lower bound sets the vector of infidelity distribution parameters (see [Definitions](#definitions) section).

#### Usage
- `d = rt_bound(dm, proto, nshots, 'state')` returns the vector `d` of the lower bound of parameters estimation variances for the quantum state density matrix `dm` and measurements protocol described by `proto` and `nshots` (see [Data format](#format) section); the purification rank is `rank(dm)`
- `d = rt_bound(chi, proto, nshots, 'process')` returns the lower bound for the quantum process matrix `chi`
- `d = rt_bound(chi, proto, nshots, 'process', 'TracePreserving', false)` considers a non trace-preserving process
- `d = rt_bound( ___ , 'Rank', r)` specifies the purification rank `r`

#### Output
- `d` &ndash; column-vector of the lower bound of parameters estimation variances

### <a name="rt_lossfun">rt_lossfun</a>
Calculates the tomography loss function for a quantum state or a quantum process (see [Definitions](#definitions) section). The result is based on the calculation of the corresponding bounded variances with [rt_bound](#rt_bound).

#### Usage
- `l = rt_lossfun(dm, proto, 'state')` returns the loss function value `l` for the quantum state density matrix `dm` and measurements operators described by `proto` (see [Data format](#format) section); the purification rank is `rank(dm)`
- `l = rt_lossfun(chi, proto, 'process')` returns the loss function value for the quantum process matrix `chi`
- `l = rt_lossfun(chi, proto, 'process', 'TracePreserving', false)` considers a non trace-preserving process
- `l = rt_lossfun( ___ , 'Rank', r)` specifies the purification rank `r`

#### Output
- `l` &ndash; loss function value

## <a name="rt_gchi2pdf">rt_gchi2pdf</a>
Calculates the generalized chi-squared distribution (see [Definitions](#definitions) section) probability density function.

#### Usage
- `p = rt_gchi2pdf(x, d)` returns the generalized chi-squared probability density `p` function with the variance vector `d` for points `x`
- `[p, x] = rt_gchi2pdf([], d)` automatically calculates the points grid and returns it as a second argument

#### Output
- `p` &ndash; probability density function values
- `x` &ndash; points of function calculation

## <a name="rt_gchi2cdf">rt_gchi2cdf</a>
Calculates the generalized chi-squared distribution (see [Definitions](#definitions) section) cumulative distribution function.

#### Usage
- `f = rt_gchi2cdf(x, d)` returns the generalized chi-squared cumulative distribution function `f` with the variance vector `d` for points `x`
- `[f, x] = rt_gchi2cdf([], d)` automatically calculates the points grid and returns it as a second argument

#### Output
- `f` &ndash; cumulative distribution function values
- `x` &ndash; points of function calculation

## <a name="rt_gchi2inv">rt_gchi2inv</a>
Calculates the generalized chi-squared distribution (see [Definitions](#definitions) section) inverse cumulative distribution function.

#### Usage
- `x = rt_gchi2pdf(f, d)` returns the generalized chi-squared inverse cumulative distribution function `x` with the variance vector `d` for points `f`

#### Output
- `x` &ndash; inverse cumulative distribution function values

### <a name="rt_randstate">rt_randstate</a>
Generates a fixed rank quantum state using the partial tracing. A Haar random pure state in the <img alt="rd" src="https://render.githubusercontent.com/render/math?math=rd" />-dimensional Hilbert space is generated and then the partial trace over <img alt="r" src="https://render.githubusercontent.com/render/math?math=r" />-dimensional subsystem is taken.

#### Usage
- `dm = rt_randstate(dim)` returns a random full-rank density matrix `dm` for the system of dimension `dim`
- `dm = rt_randstate(dim, 'Rank', r)` specifies the state rank `r`

#### Output
- `dm` &ndash; density matrix

### <a name="rt_randunitary">rt_randunitary</a>
Generates a Haar random unitary matrix.

#### Usage
- `u = rt_randunitary(dim)` returns a random unitary matrix for the system of dimension `dim`

#### Output
- `u` &ndash; unitary matrix

### <a name="rt_randprocess">rt_randprocess</a>
Generates a fixed rank quantum process using the extended dynamics representation. A Haar random unitary matrix in the <img alt="rd" src="https://render.githubusercontent.com/render/math?math=rd" />-dimensional Hilbert space is generated and then the partial trace over <img alt="r" src="https://render.githubusercontent.com/render/math?math=r" />-dimensional subsystem dynamics is taken. The process representations are described in the [Definitions](#definitions) section.

#### Usage
- `chi = rt_randprocess(dim)` returns a random rank-1 process matrix `chi` for the system of dimension `dim`
- `chi = rt_randprocess(dim, 'Rank', r)` specifies the process rank `r`
- `kr = rt_randprocess( ___ , 'Form', 'kraus')` returns the Kraus representation `kr` of the process: `kr(:,:,k)` is the k-th Kraus operator
- `__ = rt_randprocess( ___ , 'TracePreserving', false)` considers a non trace-preserving process

#### Output
- `chi` or `kr` &ndash; process matrix or Kraus operators

### <a name="rt_process_reform">rt_process_reform</a>
Changes the quantum process representation (see [Definitions](#definitions) section and [rt_randprocess](#rt_randprocess) description).

#### Usage
- `kr = rt_process_reform(chi, 'chi2kraus')` returns the process Kraus operators matrices `kr` by process matrix `chi`
- `kr = rt_process_reform(chi, 'chi2kraus', r)` specifies the number of Kraus operators `r`
- `e = rt_process_reform(chi, 'chi2root')` returns the process matrix `chi` square root
- `e = rt_process_reform(chi, 'chi2root', r)` specifies the purification rank `r`
- `chi = rt_process_reform(e, 'root2chi')` returns the process matrix `chi` by its square root `e`
- `kr = rt_process_reform(e, 'root2kraus')` returns the process Kraus operators matrices `kr` by the process matrix square root `e`
- `chi = rt_process_reform(kr, 'kraus2chi')` returns the process matrix `chi` by the corresponding Kraus operators matrices `kr`
- `e = rt_process_reform(kr, 'kraus2root')` returns the process matrix square root by the corresponding Kraus operators matrices `kr`

#### Output
- `chi`, or `kr`, or `e` &ndash; process matrix, or Kraus operators, or process matrix square root

### <a name="rt_proto_measurement">rt_proto_measurement</a>
Generates the measurements operators of a specific type in agreement with the [data format](#format).

#### Usage
- `proto = rt_proto_measurement('mub', 'Dim', dim)` returns a mutually-unbiased bases protocol `proto` for the system of dimension `dim` (only 2, 3, 4 and 8 are currently supported)
- `proto = rt_proto_measurement('tetra')` returns the qubit projective measurements protocol according to the tetrahedron symmetry in the Bloch sphere
- `proto = rt_proto_measurement('tetra', 'Modifier', 'operator+')` only +1 eigenvectors are considered
- `proto = rt_proto_measurement('tetra', 'Modifier', 'operator-')` only -1 eigenvectors are considered
- `proto = rt_proto_measurement( ___ , 'Modifier', 'operator')` each protocol operator is measured separately
- `proto = rt_proto_measurement( ___ , 'NSub', n)` considers an `n`-th tensor power of the protocol describing factorized measurements of `n` subsystems

#### Output
- `proto` &ndash; a cell array of measurements operators matrices

### <a name="rt_proto_preparation">rt_proto_preparation</a>
Generates a set of density matrices according to a specific protocol.

#### Usage
- `proto = rt_proto_preparation('mub', 'Dim', dim)` returns a set of `dim * (dim + 1)` density matrices `proto` of a mutually-unbiased bases for the system of dimension `dim` (only 2, 3, 4 and 8 are currently supported)
- `proto = rt_proto_preparation('tetra')` returns a set of `4` qubit density matrices according to the tetrahedron symmetry in the Bloch sphere
- `proto = rt_proto_preparation('octa')` returns a set of `8` qubit density matrices according to the octahedron symmetry in the Bloch sphere
- `proto = rt_proto_preparation( ___ , 'NSub', n)` considers an `n`-th tensor power of the protocol describing factorized preparation of `n` subsystems

#### Output
- `proto` &ndash; a 3D array of density matrices matrices along the third dimension

### <a name="rt_proto_process">rt_proto_process</a>
Generates a quantum process tomography protocol as a tensor product of the preparation and measurement protocols in agreement with the [data format](#format).

#### Usage
- `proto = rt_proto_process(proto_prep, proto_meas)` returns the quantum process tomography protocol `proto` as a `2`-by-`mp*mm` cell array where the first row sets the input input density matrices and the second row sets the output measurements operators; `proto_prep` is a size-`mp` preparation protocol and `proto_meas` is a size-`mm` measurement protocol (see [Data format](#format) section)

#### Output
- `proto` &ndash; a cell array with the input density matrices in the first row and the measurements operators in the second one; the columns are for different measurements schemes

### <a name="rt_nshots_divide">rt_nshots_divide</a>
Divides a total integer sample size equally over a set of measurements.

#### Usage
- `nshots = rt_nshots_divide(n, m)` returns a length `m` column-vector `nshots` of equal integers that sum up to `n`. If `n` if not divisible by `m` the first elements of `nshots` are one unit greater than the rest ones
- `nshots = rt_nshots_divide(n, m, 'total')` - makes all the elements of `nshots` strictly equal to each other by allowing non integer values
- `nshots = rt_nshots_divide(n, m, 'equal')` - makes all the elements of `nshots` equal `n` so the `nshots` sums up to `n*m`

#### Output
- `nshots` &ndash; a column-vector of measurements repetitions

### <a name="rt_meas_matrix">rt_meas_matrix</a>
Generates the measurement matrix <img alt="B" src="https://render.githubusercontent.com/render/math?math=B" /> from the measurements operators <img alt="P_j" src="https://render.githubusercontent.com/render/math?math=P_j" /> such that the _j_-th element of the column vector <img alt="B*\vec\rho" src="https://render.githubusercontent.com/render/math?math=B\cdot\vec\rho" /> (<img alt="\vec\rho" src="https://render.githubusercontent.com/render/math?math=\vec\rho" /> is the reshaped into column density matrix <img alt="rho" src="https://render.githubusercontent.com/render/math?math=\rho" />) is equal to <img alt="trace(rho*P_k)" src="https://render.githubusercontent.com/render/math?math=\textrm{Tr}(\rho P_k)" />.

#### Usage
- `B = rt_meas_matrix(P)` returns the measurement matrix `B` by the 3D array with measurements operators along the third dimension

#### Output
- `B` &ndash; measurement matrix

### <a name="rt_kron3d">rt_kron3d</a>
Calculates the generalized Kronecker product of two 3D arrays.

#### Usage
- `A = rt_kron3d(A1, A2)` returns the generalized Kronecker product `A` of 3D arrays `A1` and `A2` such that the result contains elements `kron(A1(:,:,j1), A2(:,:,j2))` with all possible combinations of `j1` and `j2` indices

#### Output
- `A` &ndash; a 3D array

### <a name="rt_pinv">rt_pinv</a>
Performs the quantum state reconstruction using Moore-Penrose pseudo-inversion.

#### Usage
- `dm = rt_pinv(P, freq)` returns a full-rank density matrix `dm` by the 3D array `P` with measurements operators along the third dimension and the column-vector `freq` with corresponding observed frequencies
- `dm = rt_pinv(P, freq, r)` specifies the density matrix rank `r`
- `[dm, c] = rt_pinv( ___ )` also returns the state square root `c`

#### Output
- `dm` &ndash; density matrix
- `c` &ndash; density matrix square root

### <a name="rt_purify">rt_purify</a>
Performs the density matrix purification by taking its square root. If the input matrix has negative eigenvalues it is projected to the set of valid density matrices.

#### Usage
- `c = rt_purify(dm)` returns the density matrix `dm` square root `c` as a matrix with `rank(dm)` columns
- `c = rt_purify(dm, r)` specifies purification rank `r`

#### Output
- `c` &ndash; density matrix square root

### <a name="rt_prttrace">rt_prttrace</a>
Calculates the partial trace of a two-component system density matrix.

#### Usage
- `dms = rt_prttrace(dm, dims, sind)` returns the partial trace `dms` of the density matrix `dm` over a subsystem with index `sind` (`1` or `2`); the subsystems dimensions array `dims = [dim1, dim2]` should have the product equal to `size(dm, 1)`

#### Output
- `dms` &ndash; subsystem density matrix

### <a name="rt_to_simplex">rt_to_simplex</a>
Calculates the projection of a real-valued vector onto standard simplex.

#### Usage
- `pp = rt_to_simplex(p)` returns the projection `pp` of the vector-column `p` onto standard simplex

#### Output
- `p` &ndash; non negative column-vector

### <a name="rt_pval">rt_pval</a>
Calculates the chi-squared test p-value for a specific number of degrees of freedom. If the case of high chi-squared values the function uses the variable-precision arithmetic.

#### Usage
- `pval = rt_pval(chi2, df)` returns the p-value of the chi-squared test for the chi-squared value `chi2` and the number of degrees of freedom `df`

#### Output
- `pval` &ndash; p-value of the chi-squared test

## <a name="rt_optimizer">rt_optimizer</a>
The class for solving optimization tasks (see [Algorithms](#algorithms) section).

### Creating an object
- `optim = rt_optimizer(type)` creates an object for a specific optimization task `type`; possible values are `'fixed_point'` for the fixed-point iteration method, `'proximal_ascend'` for the proximal gradient ascend, `'auto_rank'` for the model rank optimization

### Fields
- `fixed_point_opts` &ndash; fixed-point iterations method options
	* `fixed_point_opts.display` &ndash; display output (default: `true`)
	* `fixed_point_opts.max_iter` &ndash; maximal number of iterations (default: `1e6`)
	* `fixed_point_opts.tol` &ndash; optimization tolerance (default: `1e-8`)
	* `fixed_point_opts.reg_coeff` &ndash; regularization coefficient (default: `0.5`)
- `proximal_ascend_opts` &ndash; proximal gradient ascend method options
	* `proximal_ascend_opts.display` &ndash; display output (default: `true`)
	* `proximal_ascend_opts.max_iter` &ndash; maximal number of iterations (default: `1e6`)
	* `proximal_ascend_opts.tol` &ndash; optimization tolerance (default: `1e-8`)
- `auto_rank_opts` &ndash; model rank optimization options
	* `auto_rank_opts.display` &ndash; display output (default: `true`)
	* `auto_rank_opts.sl` &ndash; chi-squared test significance level (default: `0.05`)

### Methods
#### set_options
Sets the optimization options in the form of Name-Value pairs, where Name specifies the object field codename and the Value specifies the field value.

- `optim.set_options(code1, value1, code2, value2, ___ )`

#### run
Runs the optimization procedure with specific arguments.

- `x = optim.run( ___ )` returns the optimization result `x`; the arguments list and output depends on the optimization task (see below)
- `[x, info] = optim.run( ___ )` also returns addition optimization information `info`

### Fixed-point iteration method
#### Arguments
- `x0` &ndash; starting point (a general complex-valued ND array)
- `fVal` &ndash; handler of fixed-point function value (the equation to solve is `fVal(x) = x`)

#### Output
- `x` &ndash; optimization result
- `info.iter` &ndash; number of iterations

### Proximal gradient ascend method
#### Arguments
- `x0` &ndash; starting point (a general complex-valued ND array)
- `fFun` &ndash; handler for function value computation (the function value for the point `x` is `fFun(x)`)
- `fdFun` &ndash; handler for function gradient computation (the function gradient for the point `x` is `fdFun(x)`)
- `fProx` &ndash; handler for proximal operator (the result of proximal operator action on the point `x` is `fProx(x)`)
- `LipsConst` &ndash; Lipschitz constant of the function

#### Output
- `x` &ndash; optimization result
- `info.iter` &ndash; number of iterations

### Model rank optimizer
#### Arguments
- `rmax` &ndash; maximal rank value
- `fData` &ndash; handler for calculating data: `fData(r)` should be a structure array containing p-value field `pval` of rank-`r` model

#### Output
- `x` &ndash; data for optimal rank model
- `info{r}` &ndash; data for rank-`r` model

## <a name="references">References</a>
<a name="ref1">[1]</a> Bogdanov Yu. I. Unified statistical method for reconstructing quantum states by purification // _JETP **108(6)**_, 928-935 (2009); doi: 10.1134/S106377610906003X

<a name="ref2">[2]</a> Huan L. and Zhouchen L. Accelerated proximal gradient methods for nonconvex programming // Advances in Neural Information Processing Systems 28 (NIPS 2015)
