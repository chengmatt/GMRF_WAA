## Unlocking the triad of age, year, and cohort effects in stock assessment: a proof-of-concept study

### Authors
Matthew LH. Cheng <sup>1</sup>

<sup>1</sup> Department of Fisheries at Lena Point, College of Fisheries and Ocean Sciences, University of Alaska Fairbanks, 17101 Point Lena Loop Rd, Juneau, Alaska 99801, USA

### Abstract
Many demographic processes vary by age and over time, and accounting for this variation within fisheries management remains a key challenge for many contemporary stock assessments. Although there is evidence for time, age, and cohort specific effects on various components within stock assessment (e.g., selectivity, growth), methods are lacking to simultaneously estimate autocorrelation over time, among ages, and by cohort while also quantifying residual variation. Drawing from previous research on separable cohort models, we reintroduce the idea of “triple-separability”, which simultaneously estimates autocorrelation for time, age, and cohort effects, and reduces to two-dimensional autocorrelation when one of the processes is fixed at zero. Utilizing eastern Bering Sea walleye pollock (Gadus chalcogrammus) as a case-study, we illustrate differences in predicted weight-at-age values from models with and without a triple-separable assumption and show that traditional model selection tools can be used to identify the relative evidence and magnitude of age, time, and cohort effects. We recommend that the method be integrated as a routine option within next-generation stock assessments, and note that it generalizes widely-used options in existing assessment frameworks (i.e., WHAM and SAM).
### Keywords
mixed-effects, stock assessment, state-space, age-and time-varying, cohort effects

### Repository Structure
| Folder  | Items |
| --------| --------|
|data| Weight-at-age matrix with standard deviations for EBS walleye pollock |
|docs| Auxiallary documents describing triple-separability in the VPA context |
|figs| Contains general figures and model outputs |
|R_scripts| Contains scripts for model runs, visualizations, and constructing a precision matrix with marginal variance |
|src| Contains source code for TMB model, which is compiled through R|
|test|R scripts testing the development of the TMB model and checks for correct model dimensions|
|output| Model outputs from 2x2x2 factorial experiment |

### Specification of the Precision Matrix
| Purpose  | Location |
| --------| --------|
|Precision Matrix| Construction of the precision matrix with marginal variance can be [found here](https://github.com/chengmatt/Triple_Separability/blob/master/R_scripts/make_precision/Construct_precision_2023-01-02.R) as a R function, or [here](https://github.com/chengmatt/Triple_Separability/blob/fa8e8dfe0e44b29f5fd98726352fe50aea8e9db8/src/triple_sep_waa.cpp#L8-L126?plain=1) as a local function defined in TMB|
|Condtional Variance &Omega; matrix | Code for implementing the &Omega; matrix with condtional variance can be found [here](https://github.com/chengmatt/Triple_Separability/blob/6ebbeecfe4dd002de66c655d7d687f688cdd1954/src/triple_sep_waa.cpp#L72-L81?plain=1)|
|Marginal Variance &Omega; matrix | Code for implementing the &Omega; matrix with marginal variance in TMB can be found [here](https://github.com/chengmatt/Triple_Separability/blob/6ebbeecfe4dd002de66c655d7d687f688cdd1954/src/triple_sep_waa.cpp#L83-L116?plain=1)|

