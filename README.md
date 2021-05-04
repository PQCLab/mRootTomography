# Root approach quantum tomography

MATLAB library for the discrete variables quantum state and quantum process tomography using root approach. The library contains a set of tools for quantum state and quantum process reconstruction by the complementary measurements results, estimation of statistical adequacy and theoretical analysis of reconstruction fidelity. Full product documentation is available [here](Documentation.md).

- [Getting Started](#start)
- [Definitions, algorithms and data format](#format)
- [Quantum state tomography](#qst)
- [Quantum process tomography](#qpt)
- [License](#license)

## <a name="start">Getting Started</a>

Required toolboxes:
* Statistics toolbox
* Symbolic Math Toolbox (optional)

To install the library clone the repository or download and unpack zip-archive. Run the startup script before use.

```
>> rt_startup
```

## <a name="format">Definitions, algorithms and data format</a>

For the entities definitions, algorithms description and required data format see full product [documentation](Documentation.md).

## <a name="qst">Quantum state tomography</a>

Examples directory of the project contains a set of examples that show basic features of the library. Below we briefly review the quantum state tomography example.

Consider an example of a pure ququart state reconstruction using mutually-unbiased bases measurement protocol.
```
dim = 4; % System dimension
r_true = 1; % True state rank
dm_true = rt_randstate(dim, 'Rank', r_true); % True state
nshots = 1e3; % Total sample size
proto = rt_proto_measurement('mub', 'dim', dim); % Generate measurements operators
```

The `rt_experiment` class allows one to store and simulate tomography data.
```
ex = rt_experiment(dim, 'state'); % Generate experiment instance
ex.set_data('proto', proto, 'nshots', nshots); % Set measurements protocol
ex.simulate(dm_true); % Simulate measurement data
```

The reconstruction is performed using the `rt_dm_reconstruct` function. By default the state rank is estimated automatically using the adequacy criteria.
```
dm_rec = rt_dm_reconstruct(ex, 'Display', 10);
Fidelity = rt_fidelity(dm_rec, dm_true);
fprintf('Fidelity: %.6f\n', Fidelity);
```

Output:
```
=== Automatic rank estimation ===
=> Try rank 1
Optimization: fixed point iteration method
Iteration 33 		 Delta 8.0049e-09
=> Rank 1 is statistically significant at significance level 0.05. Procedure terminated.
Fidelity: 0.996887
```

Using the fiducial approach and the theoretical infidelity distribution one can use the reconstruction result to estimate the guaranteed reconstruction fidelity. In the following code we estimate the 95%-probability fidelity bound ![F_95](https://latex.codecogs.com/svg.latex?F_%7B95%7D). That means that we get the fidelity ![F_95](https://latex.codecogs.com/svg.latex?F_%7B95%7D) or higher with probability 95%.
```
d = rt_bound(dm_rec, ex); % Calculate variances
Fidelity95 = 1 - rt_gchi2inv(0.95, d); % Get fidelity bound
fprintf('Fiducial 95%% fidelity bound: %.6f\n', Fidelity95);
```

Output:
```
Fiducial 95% fidelity bound: 0.992890
```

The following code plots the infidelity distribution based on the true state and shows the fidelity of reconstructed state.
```
d = rt_bound(dm_true, ex);
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');
```

![Theoretical distribution and reconstruction result](Examples/infiddistr.png?x=2)

The theoretical infidelity mean and variance could be obtained from the variance vector.
```
sum(d) % Mean infidelity
2*sum(d.^2) % Infidelity variance
1-sum(d) % Mean fidelity
```

## <a name="qpt">Quantum process tomography</a>
The measurements protocol for the quantum process tomography could be generated by the combination of preparation and measurements protocols.
```
proto_prep = rt_proto_preparation('tetra');
proto_meas = rt_proto_measurement('mub', 'dim', dim);
proto = rt_proto_process(proto_prep, proto_meas);
```

The reconstruction of a quantum process is performed in a very similar way to the quantum state reconstruction.
```
chi = rt_chi_reconstruct(ex);
```

The fidelity bound estimation is also done in a similar way.
```
d = rt_bound(chi, ex);
```

## <a name="license">License</a>

All code found in this repository is licensed under GPL v3
