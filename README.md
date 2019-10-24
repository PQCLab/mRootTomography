# Root approach quantum tomography

MATLAB library for the discrete variables quantum state tomography using root approach. The library contains a set of tools for quantum state reconstruction by the complementary measurements results, estimation of statistical adequacy and theoretical analysis of reconstruction fidelity. Full product documentation is available [here](Documentation.md).

- [Getting Started](#start)
	* [Prerequisites and installing](#install)
	* [Data format](#format)
	* [Quantum state tomography](#qst)
	* [Fidelity distribution](#fidelity)
- [Algorithms](#algorithms)
- [References](#references)

## <a name="start">Getting Started</a>

<a name="install" />
### Prerequisites and installing

Required toolboxes:
* Statistics toolbox

To install the library clone the repository or download and unpack zip-archive. Before using run the startup script.

```
>> rt_startup
```

<a name="format" />
### Data format

Measurement protocol is defined by a set of complementary measurement experiments over density matrix ![rho](https://latex.codecogs.com/svg.latex?%5Crho). Every experiment is repeated many times and has several possible measurement outcomes. The probability to get ![k](https://latex.codecogs.com/svg.latex?k)-th outcome is determined by the measurement operator ![M_k](https://latex.codecogs.com/svg.latex?M_k) as ![p_k=trace(rho*M_k)](https://latex.codecogs.com/svg.latex?p_k%3D%5Ctext%7BTr%7D%28%5Crho%20M_k%29). The set of measurement operators and the number of experiments repetitions define the **_measurement protocol_**. The number of observations for each outcome define the **_measurement results_**. The following code describe the required data format.
```
proto{j}(:,:,k) % Measurement operator matrix for k-th outcome in j-th experiment
nshots(j) % Number of j-th experiment repetitions
clicks{j}(k) % Number of k-th outcome observations in j-th experiment
```

The following code generates an example of ideal ![N](https://latex.codecogs.com/svg.latex?N)-qubit Pauli measurements protocol (each qubit is measured in ![sigma_x](https://latex.codecogs.com/svg.latex?%5Csigma_x), ![sigma_y](https://latex.codecogs.com/svg.latex?%5Csigma_y) and ![sigma_z](https://latex.codecogs.com/svg.latex?%5Csigma_z) bases) where every of ![3^N](https://latex.codecogs.com/svg.latex?3%5EN) measurement experiments is repeated ![n](https://latex.codecogs.com/svg.latex?n) times.
```
N = 1; % Number of qubits
n = 1e3; % Number of repetitions in every experiment
proto = rt_proto_measurement('pauli', N);
nshots = rt_nshots_devide(n,length(proto),'equal');
```

You can simulate the measurements of some density matrix.
```
dm_true = rt_randstate(2,'mixed',1); % Generate random one-qubit pure state density matrix
clicks = rt_simulate(dm_true, proto, nshots); % Simulate measurements
```

<a name="qst" />
### Quantum state tomography
Reconstruct density matrix and estimate fidelity comparing to expected state.
```
dm_expected = [0.5, 0.45; 0.45, 0.5]; % Expected density matrix
dm = rt_dm_reconstruct(clicks,proto,nshots,'Display',true);
fprintf('Fidelity: %.6f\n', rt_fidelity(dm, dm_expected));
```
Output:
```
=== Automatic rank estimation ===
Try rank 1
Iteration 311 		 Difference 9.8078e-09
Try rank 2
Iteration 38 		 Difference 8.4052e-09
Rank 2 is statistically significant at significance level 0.005. Procedure terminated.
Fidelity: 0.999641
```
If you don't specify `nshots` then number of ![j](https://latex.codecogs.com/svg.latex?j)-th experiment runs is taken as the sum over all counts.
```
rt_dm_reconstruct(clicks,proto) % nshots(j) = sum(clicks{j})
```
Instead of automatic rank estimation you can pick a specific rank of the quantum state model. For example, in some experiments it could be likely to have a pure (rank 1) state.
```
dm_expected = [1, 1j; -1j, 1]/2;
dm = rt_dm_reconstruct(clicks,proto,nshots,'Display',true,'Rank',1);
fprintf('Fidelity: %.6f\n', rt_fidelity(dm, dm_expected));
```
Output:
```
Iteration 20 		 Difference 8.5074e-09
Fidelity: 0.998998
```
In general, you can specify any rank from 1 to the Hilbert space dimension.

<a name="fidelity" />
### Fidelity distribution

One can theoretically estimate reconstruction fidelity for the desired measurement and quantum state model. This code will return vector of variances of the quantum state independent parameters.
```
d = rt_dm_theory(dm_expected,proto,nshots)
```
The distribution of infidelity has the form of a generalized chi-squared distribution: ![1-F](https://latex.codecogs.com/svg.latex?%5Cinline%201-F%20%5Csim%20%5Csum_j%7Bd_j%20%5Cxi_j%5E2%7D),  where each ![xi_j](https://latex.codecogs.com/svg.latex?%5Cxi_j) is independent and has a normal distribution with zero mean and unit variance [[1]](#ref1). The following script plots theoretical distribution and shows the value of infidelity for a single tomography experiment.

```
dm_expected = [0.5, 0.45; 0.45, 0.5];
dm = rt_dm_reconstruct(clicks,proto,nshots); % Reconstruct state
F = rt_fidelity(dm, dm_expected); % Estimate fidelity with expected state
d = rt_dm_theory(dm_expected,proto,nshots); % Find theoretical distribution parameters
[p,dF] = rt_chi2pdf_general([],d); % Get distribution
figure; hold on; grid on;
plot(dF,p,'LineWidth',1.5,'DisplayName','Theory');
plot([1-F,1-F],ylim,'LineWidth',1.5,'DisplayName','Reconstruction');
xlabel('$$1-F$$','Interpreter','latex');
legend('show');
```

![Theoretical distribution and reconstruction result](Examples/infiddistr.png)

The theoretical infidelity mean and variance could be obtained from the variance vector.
```
sum(d) % Mean infidelity
2*sum(d.^2) % Infidelity variance
1-sum(d) % Mean fidelity
```
<a name="algorithms" />
## Algorithms

Consider a quantum state in the Hilbert space ![H](https://latex.codecogs.com/svg.latex?%5Cmathcal%7BH%7D) described by the density matrix ![rho](https://latex.codecogs.com/svg.latex?%5Crho).
According to the root approach the rank-![r](https://latex.codecogs.com/svg.latex?r) quantum state parametrization is defined by the complex matrix ![c](https://latex.codecogs.com/svg.latex?c) of size ![sxr](https://latex.codecogs.com/svg.latex?s%20%5Ctimes%20r) (![s=dimH](https://latex.codecogs.com/svg.latex?s%3D%5Ctext%7Bdim%7D%5Cmathcal%7BH%7D)) such that ![rho=cc+](https://latex.codecogs.com/svg.latex?%5Crho%3Dcc%5E%5Cdagger). To get the quantum state maximum likelihood estimation one must solve the following quasi-linear equation (_likelihood equation_) [[1]](#ref1):
<p align="center"><img src="https://latex.codecogs.com/svg.latex?Ic%3DJ%28c%29c"/></p>
where
<p align="center"><img src="https://latex.codecogs.com/svg.latex?I%3D%5Csum_%7Bj%7D%7Bn_j%20M_j%7D%2C%5C%3B%5C%3B%5C%3B%5C%3BJ%28c%29%3D%5Csum_%7Bj%7D%7B%5Cfrac%7Bk_j%7D%7B%5Ctext%7BTr%7D%28cc%5E%5Cdagger%20M_j%29%7D%20M_j%7D"/></p>

The sums here are taken over all measurement experiments and possible outcomes in them. ![k_j](https://latex.codecogs.com/svg.latex?k_j) is the number of observed outcomes corresponding to the measurement operator ![M_j](https://latex.codecogs.com/svg.latex?M_j) and the number of measurements repetitions ![n_j](https://latex.codecogs.com/svg.latex?n_j).

The search of the likelihood equation solution is performed by the fixed-point iteration method:
<p align="center"><img src="https://latex.codecogs.com/svg.latex?c_%7Bi&plus;1%7D%3D%281-%5Calpha%29I%5E%7B-1%7DJ%28c_i%29c_i%20&plus;%20%5Calpha%20c_i"/></p>

Here ![alpha](https://latex.codecogs.com/svg.latex?%5Calpha) is the regularization parameter. We use Moore-Penrose pseudo-inversion to get ![c_0](https://latex.codecogs.com/svg.latex?c_0).

<a name="references" />
## References

<a name="ref1">[1]</a> Bogdanov Yu. I. Unified statistical method for reconstructing quantum states by purification // _JETP_ **108(6)** 928-935 (2009); doi: 10.1134/S106377610906003X