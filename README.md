# U1_model
One plaquette U(1) model through Lefschetz thimble method - Reproduction of arXiv:1308.0233

This python3 code reproduce the results obtain in arXiv:1308.0233.

The main part is the U1.py. It reads the parameters from param.py.
To be launch, it requires two arguments which are N_tau and beta (see the paper for details)

python3 U1.py Nt beta

The code produces gauge field configurations and measurement of the phase e^{i\phi}.

The data are analysed through the binning method with error_computation.py, through Jackknife
and Gamma method with jackknife.py.

The file plots.py produces results graph.
