# CoMMSur


**CoMMSur** implements **CoMM-Surv**, a Bayesian survival mediation framework for identifying mediators between an exposure and a time-to-event outcome.

The method jointly models exposure–mediator and mediator–outcome pathways using a **configuration-based prior**. By placing a structured prior on the pair of coefficients corresponding to the exposure–mediator and mediator–outcome effects, CoMMSur estimates the posterior probability that each mediator is actively involved in the mediation pathway.


The model is implemented in **Stan** .

---

## Installation

```r
install.packages("remotes")
remotes::install_github("yhc0211/CoMMSur")