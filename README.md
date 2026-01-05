# Interference with Systemic Negative Feedback as a Potential Mechanism for Nonmonotonic Dose-Responses of Endocrine-Disrupting Chemicals

Zhenzhen Shi<sup>1,3</sup>, Shuo Xiao<sup>2</sup>, and Qiang Zhang<sup>1</sup> 

Published in Toxicological Sciences (2025) https://doi.org/10.1093/toxsci/kfaf060

1. Gangarosa Department of Environmental Health, Rollins School of Public Health, Emory University, Atlanta, GA 30322, USA

2. Department of Pharmacology and Toxicology, Ernest Mario School of Pharmacy, Environmental and Occupational Health Sciences Institute (EOHSI), Center for Environmental Exposures and Disease (CEED), Rutgers University, Piscataway, NJ 08854, USA

3. Present address: Stritch School of Medicine, Loyola University Chicago, Chicago, IL, USA.
 

**Abstract:**
Environmental endocrine-disrupting chemicals (EDCs) often exhibit nonmonotonic dose-response (NMDR) relationships, posing significant challenges to health risk assessment and regulations. Several molecular mechanisms operating locally in cells have been proposed, however, whether and how systemic negative feedback – a global structure of all homeostatic endocrine systems – may render NMDRs is poorly understood. We hypothesized that an EDC may produce nonmonotonic effects by competing with the endogenous hormone for receptors simultaneously (i) at the central site to interfere with the feedback regulation, and (ii) at the peripheral site to disrupt the hormone’s endocrine action. We constructed a dynamical model of a generic hypothalamic-pituitary-endocrine (HPE) axis with negative feedback to evaluate the hypothesis and biological conditions that favor NMDR. Our modeling found that when an EDC interferes sufficiently with the central feedback action, the net endocrine effect at the peripheral target site can be opposite to what is expected of an agonist or antagonist at low concentrations. J/U or Bell-shaped NMDRs arise when the EDC has differential binding affinities and/or efficacies, relative to the endogenous hormone, for the peripheral and central receptors. Novel quantitative relationships between these biological parameter variabilities and associated distributions were discovered, which can distinguish J/U and Bell-shaped NMDRs from monotonic responses. In conclusion, the ubiquitous negative feedback regulation in endocrine systems may act as a universal mechanism for counterintuitive and nonmonotonic effects of EDCs. Depending on key receptor kinetic and signaling properties of EDCs and endogenous hormones, certain individuals may be more susceptible to these complex endocrine effects.

**Keywords:** endocrine-disrupting chemicals, nonmonotonic dose-response, negative feedback, binding affinity, efficacy


#  MATLAB Model Code and Simulated Datasets
- EDC_cmd.m: Main MATLAB file simulating the model and producing figures
- EDC_ode.m: ODE file of the HPE model
- Default_param.mat: Default model parameter values
- 6_parameter_MC_simulation_results/: Contains the file for 6-parameter MC simulation results
- Population_MC_simulation_results/: Contains the files for population MC simulation results
