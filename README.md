# fourierPopulationInference
Infer Fourier-based population function from phylogenetic trees

- FourierPopulation.java and TruncatedFourierSeries.java should be saved under /beast2/src/beast/base/evolution/tree/coalescent.
- All XML's use Remaster to simulate trees and incorporate with BEAST to do the inference.
    - constPop.xml, expPop.xml and periodicPop.xml use Remaster to simulate different trees based on constant, exponential, periodic population functions.
    - sampleDiffTime.xml simulates a tree in which leaves are sampled at different times.
    - fourierPop_fixTree.xml infers the Fourier population function from fixed trees.
    - BSP_fixTree.xml infers the Bayesian Skyline Plot from fixed trees.
