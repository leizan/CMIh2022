The code provided in this package supports the research paper below.

Methods:
- mixedCmiIEstimator.R contains the code for the hybrid estimator;  
- method.R contains the code for the adaptive independent test and other tests compared in experiments;
- MS-R.R contains the code for k-nearest-neighbour based CMI estimator;   
- folder 'algorithm' contains the code for histogram based CMI estimator.

Experiments:
The 'experiments' folder contains:
- syntheticGenerator.R generates ground-truth data for CMI estimators;
- test_distrib.R evaluates CMI estimators on ground-truth data (Figure 1);
- test_distrib_multidim.R evaluates CMI estimators on ground-truth data (Figure 2);
- data.R generates synthetic data for the conditional independent test;
- testSyn.R evaluates different test methods on synthetic data(Table 1);
- testReal.R evaluates different test methods on real data(Table 2).

To run the experiments, stay in this folder CMIh2022/ with your console and run, e.g.,
Rscript experiments/test_distrib.R
The real data for the conditional independent test is in folder 'ordinal_data'. The results will be stored into the 'results' folder and the p-Value for the independent test on real data will be stored into the 'pValue' folder. You may also need to install additional R packages (see requirements in 'source.R') 


