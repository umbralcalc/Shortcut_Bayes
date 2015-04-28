#### Purpose-built research code
Please feel free to contact me if you have any suggestions or questions at umbralcalc@gmail.com. To find out more about my research, please visit the link [here](https://sites.google.com/site/umbralcalc/home). The code here is all authored by myself and is therefore not especially _flashy_, however it all should perform relatively bug-free. You can find in this repository:

1. _Shortcut_Bayes.cpp_: This takes as input the Markov-Chain Monte-Carlo posterior samples which you can find [here](http://pla.esac.esa.int/pla/) obtained using CosmoMC (a Fortran Metropolis-Hastings algorithm used prolifically in Cosmology, see [here](http://cosmologist.info/cosmomc)) by the _Planck Collaboration_ and outputs approximate Bayesian evidence ratios and maximum likelihood estimates. The approximation is made by obtaining a kernel density smoothed likelihood function over 3 core observational features of the Cosmic Microwave Background: the scalar power spectrum amplitude, the scalar spectral index and the tensor-to-scalar ratio while also including the degree of Non-Gaussinity. The code takes its priors from a class of multi-field slow-roll cosmic inflation models and compares them to predictions of some single-field mononomial potentials.
