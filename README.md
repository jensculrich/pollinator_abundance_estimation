# mark_recapture
## estimation of population size for bees in urban parks using mark-recapture detection history data
## started May 25, 2022
#### the model builds from an M0 type model (constant individual and temporal detection probabilities) described by Kery and Schaub in "Bayesian population analysis using WinBUGS: a hierarchical perspective", however, using STAN rather than WinBUGS, and will ultimately build towards considering:
#### 1) either two categorical site types or sites varying continuously in their plant community composition 
#### 2) multiple species and two years worth of data (although likely not integrated in a hierarchical way, given that the number of species considered (3 or 4) and years considered (2) are low.
