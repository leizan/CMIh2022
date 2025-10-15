# Code Python/R CMIh and the corresponding independent test

The code provided in this package supports the research paper below:
Zan, Lei, Anouar Meynaoui, Charles K. Assaad, Emilie Devijver, and Eric Gaussier. 2022. "A Conditional Mutual Information Estimator for Mixed Data and an Associated Conditional Independence Test" Entropy 24, no. 9: 1234. https://doi.org/10.3390/e24091234

## Note:

This method has also been implemented in the **Tigramite** library — a powerful toolkit for causal discovery in time series.
Special thanks to the Tigramite team for their excellent work and contribution.
👉 [https://github.com/jakobrunge/tigramite](https://github.com/jakobrunge/tigramite)


## Methods:
- mixedCmiIEstimator.R contains the code for the hybrid estimator;  
- method.R contains the code for the adaptive independent test and other tests compared in experiments;
- MS-R.R contains the code for k-nearest-neighbour based CMI estimator;   
- folder 'algorithm' contains the code for histogram based CMI estimator.

## Experiments:
The 'experiments' folder contains:
- syntheticGenerator.R generates ground-truth data for CMI estimators;
- test_distrib.R evaluates CMI estimators on ground-truth data (Figure 1 & Figure 2(middle, right));
- test_distrib_multidim.R evaluates CMI estimators on ground-truth data (Figure 2(left) & Table 1);
- data.R generates synthetic data for the conditional independent test;
- testSyn.R evaluates different test methods on synthetic data(Table 2);
- testTemperature.R evaluates different test methods on perprocessed DWD dataset(Table 3);
- testADHD.R evaluates different test methods on ADHD-200 dataset(Table 4);
- testReal.R evaluates different test methods on EasyVista dataset(Table 5).

###
1. To run the experiments, stay in this folder CMIh2022/ with your console and run:
```bash
Rscript experiments/test_distrib.R
```
2. The real data for the conditional independent test is in folder 'ordinal_data'. The results will be stored into the 'results' folder and the p-Value for the independent test on real data will be stored into the 'pValue' folder. 

3. You may also need to install additional R packages (see requirements in 'source.R') 

4. The python code is in the folder _python_code_.
