## Estimating age, year, and cohort effects in stock assessments: demonstration of a computationally efficient and reproducible framework

### Authors
Matthew LH. Cheng,<sup>1</sup> James T. Thorson,<sup>2</sup> James N. Ianelli,<sup>3</sup> Curry J. Cunningham<sup>1</sup> 

<sup>1</sup> Department of Fisheries at Lena Point, College of Fisheries and Ocean Sciences, University of Alaska Fairbanks, 17101 Point Lena Loop Rd, Juneau, Alaska 99801, USA

<sup>2</sup>  Habitat and Ecological Processes Research Program, Alaska Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand Point Way N.E., Seattle, WA 98115, USA

<sup>3</sup>  Resource Ecology and Fisheries Management Division, Alaska Fisheries Science Center, NOAA, 7600 Sand Point Way N.E., Building 4 Seattle, WA 98115, USA



### Abstract
Many demographic processes vary by age and over time, and accounting for this variation within fisheries management remains a key challenge for contemporary stock assessments. Although there is evidence for time, age, and cohort specific effects on various components within stock assessment (e.g., selectivity, growth), methods are lacking to simultaneously estimate autocorrelation over time, among ages, and by cohort while also quantifying residual variation. We develop an approach which facilitates the simultaneous estimation of autocorrelation for time, age, and cohort effects. Utilizing eastern Bering Sea walleye pollock (Gadus chalcogrammus) as a case-study, we conduct a factorial experiment and demonstrate differences in predicted weight-at-age values from models that estimate different combinations of correlation parameters. We show that traditional model selection tools can be used to identify the relative evidence and magnitude of age, time, and cohort effects. We show that this method be easily integrated as a routine option within next-generation stock assessments, and generalizes widely-used options in existing assessment frameworks (e.g., WHAM and SAM).  Code is available in a GitHub repository that replicates our analysis and includes a Template Model Builder function for assembling the precision matrix for easy adoption in other software packages.


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
|Precision Matrix| Construction of the precision matrix with marginal variance can be [found here](https://github.com/chengmatt/GMRF_WAA/blob/master/R_scripts/make_precision/Construct_precision.R) as a R function, or [here](https://github.com/chengmatt/GMRF_WAA/blob/26b38d6e5dd64bd79d16f1c1d4f3671364bfb93d/src/GMRF_WAA.cpp#L14-L132) as a local function defined in TMB|
|Condtional Variance &Omega; matrix | Code for implementing the &Omega; matrix with condtional variance can be found [here](https://github.com/chengmatt/GMRF_WAA/blob/26b38d6e5dd64bd79d16f1c1d4f3671364bfb93d/src/GMRF_WAA.cpp#L78-L87)|
|Marginal Variance &Omega; matrix | Code for implementing the &Omega; matrix with marginal variance in TMB can be found [here](https://github.com/chengmatt/GMRF_WAA/blob/6ebbeecfe4dd002de66c655d7d687f688cdd1954/src/triple_sep_waa.cpp#L83-L116?plain=1)|

