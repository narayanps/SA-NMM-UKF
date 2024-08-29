 # Overview
 This repository contains the code and analysis used in our study on the parameter sensitivity of the Jansen and Rit Neural Mass Model (JR-NMM) for estimating neural parameters from EEG data. The JR-NMM is a well-established model for simulating the mesoscopic and macroscopic dynamics of electroencephalography (EEG), including epileptic EEG. However, the high dimensionality of the model and the uncertainty in estimating its parameters pose significant challenges. This project aims to identify the most sensitive parameters within the JR-NMM and develop a reliable method for their estimation.

 ## Repository structure and usage
 1. /func contains all the necessary functions to perform sensitivity analysis and parameter estimation
 2. run_jr_sim.m  - runs JR-NMM simulation for a given combination of input parameters
 3. run_morris_analysis - performs morris SA
 4. run_sonol_analysis - performs sobol SA
 5. test_real_data_human - estimates the parameters B and b from epileptic datasets. The datasets can be down loaded from http://ieeg-swez.ethz.ch/

SAFE toolbox (https://safetoolbox.github.io/) is needed to perform the morris analysis which we have modified to incorporate parallelization . A part of CODES toolbox (https://codes.arizona.edu/toolbox/help/html/sobol.html) has been used/modified to perform sobol analysis.
 ## Key Contributions

1. **Parameter Sensitivity Analysis**  
   We conducted a comprehensive parameter sensitivity analysis of the JR-NMM using global sensitivity analysis methods, specifically the Morris method and Sobol method. The analysis identified two key parameters, the average inhibitory synaptic gain (B) and the reciprocal of the time constant of the average inhibitory post-synaptic potentials (PSPs) (b), as the most sensitive parameters that significantly influence the model's output.

2. **Bayesian Parameter Estimation**  
   We developed a Bayesian approach for estimating the JR-NMM states and parameters based on the Expectation-Maximization (EM) algorithm combined with the Unscented Kalman Smoother (UKS-EM). This method allows for accurate estimation of the most sensitive parameters, B and b, even under varying levels of measurement noise.

3. **Application to Epileptic EEG Data**  
   The UKS-EM algorithm was applied to intracranial EEG data from 16 epileptic patients. Our results demonstrate that the parameters B and b change significantly across different seizure states, indicating their potential use in seizure tracking and control. Specifically:
   - Transition to seizure is characterized by a decrease in average B.
   - High-frequency activity during seizures is associated with an increase in b.

MIT License

Copyright (c) 2024 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

