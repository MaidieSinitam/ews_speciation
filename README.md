# Early warning signals of impending speciation

Maidie Sinitambirivoutin (1), Patrik Nosil (2), Samuel Flaxman (4), Jeffrey Feder (5), Zachariah Gompert (6), Vasilis Dakos (1,3)


1 Institute of Ecology and Environmental Sciences, Sorbonne University/CNRS/INRA/IRD/UPEC/Paris-Diderot University, Paris, France 

2 CEFE, Univ Montpellier, CNRS, EPHE, IRD, Univ Paul Valery Montpellier 3, Montpellier, 34293, France 

3 Institut des Sciences de l’Evolution Montpellier, Univ Montpellier/CNRS/EPHE/IRD, Montpellier, France 

4 Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA 5 Department of Biological Sciences, University of Notre Dame, Notre Dame, Indiana, USA

6 Department of Biology, Utah State University, Logan, UT 84322, USA 


## Abstract

Species formation is a central topic in biology and a large body of theoretical work has explored the conditions under which speciation occurs, including whether speciation dynamics are gradual or abrupt. In some cases of abrupt speciation, differentiation slowly builds up until it reaches a threshold, at which point linkage disequilibrium (LD) and divergent selection enter a positive feedback loop that triggers accelerated change. Notably, such abrupt transitions powered by positive feedbacks have also been observed in a range of other systems. Efforts to anticipate abrupt transitions have led to the development of ‘early warning signals’ (EWS), i.e. specific statistical patterns preceding abrupt transitions. Examples of EWS are rising autocorrelation and variance in time-series data, due to the reduction of the ability of the system to recover from disturbances. Here, we investigate whether speciation dynamics in theoretical models also exhibit EWS. Using a model of genetic divergence between two populations, we search for EWS before gradual and abrupt speciation events. We do so using six different metrics of differentiation: the effective migration rate, the number of selected loci, the mean fitness of our studied population, LD, FST and Dabs, a metric analogous to DXY . We find evidence for EWS, but with a heterogeneity in their strength among differentiation metrics. We specifically identify FST and the effective migration rate as the most reliable EWS of upcoming abrupt speciation events. Our results provide initial insights into potential EWS of impending speciation and contribute to efforts to generalize the mechanisms underlying EWS.

## Data

We analysed simulation outputs from the Build-up-to-Speciation ($BU2S$) model. The raw outputs are available at: <https://doi.org/10.5061/dryad.5qfttdz7d> [link](https://doi.org/10.5061/dryad.5qfttdz7d)

## Code

* `preliminary_analysis.R`: code to retrieve the relevant metrics of differentiation from the raw simulation outputs and calculate the correlation between the metrics
* `analysis_dataser.R`: code to calculate the moving autocorrelation-at-lag-1 and standard deviation of the differentiation metrics. These calculations needs a lot of computing power.
* `ews_analysis.R`: code to calculate the significance of the results using the area under the curve (AUC).

