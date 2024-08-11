# Unlocking the triad of age, year, and cohort effects for stock assessment: demonstration of a computationally efficient and reproducible framework using weight-at-age 
https://doi.org/10.1016/j.fishres.2023.106755


# Authors
Matthew LH. Cheng,<sup>1</sup> James T. Thorson,<sup>2</sup> James N. Ianelli,<sup>3</sup> Curry J. Cunningham<sup>1</sup> 

<sup>1</sup> Department of Fisheries at Lena Point, College of Fisheries and Ocean Sciences, University of Alaska Fairbanks, 17101 Point Lena Loop Rd, Juneau, Alaska 99801, USA 

<sup>2</sup>  Habitat and Ecological Processes Research Program, Alaska Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand Point Way N.E., Seattle, WA 98115, USA

<sup>3</sup>  Resource Ecology and Fisheries Management Division, Alaska Fisheries Science Center, NOAA, 7600 Sand Point Way N.E., Building 4 Seattle, WA 98115, USA

Contact: lhcheng@alaska.edu, mat.cheng412@gmail.com

# Abstract
Many demographic processes vary by age and over time, but are also hypothesized to exhibit cohort-specific patterns in variation; accounting for this variation within fisheries management remains a key challenge for contemporary stock assessments. Although there is evidence for time, age, and cohort-specific patterns in the variation of various components within stock assessment (e.g., selectivity, growth), methods are sparsely documented or are lacking to simultaneously estimate autocorrelation over time, among ages, and by cohort while also quantifying residual variation. We demonstrate an approach that facilitates the simultaneous estimation of autocorrelation for time, age, and cohort correlations, and provide two options to estimate the pointwise variance of this process (termed conditional and marginal variance). Utilizing eastern Bering Sea walleye pollock (Gadus chalcogrammus) as a case-study, we develop factorial model formulations to demonstrate differences in predicted weight-at-age values from models that estimate different combinations of correlation parameters along three axes (age, year, cohort). We show that traditional model selection tools can be used to identify the relative evidence for, and magnitude of, age, time, and cohort correlations, and demonstrate that this method can be easily integrated as a routine option within next-generation stock assessments. Code to replicate our analysis is provided in a GitHub repository (https://github.com/chengmatt/GMRF_WAA), and includes a Template Model Builder function for assembling the precision matrix to facilitate easy adoption in other software packages.


# Keywords
mixed-effects, stock assessment, state-space, age-and time-varying, cohort effects

# Repository Structure
| Folder  | Items |
| --------| --------|
|data| Weight-at-age matrix with standard deviations for EBS walleye pollock |
|docs| Auxiallary documents describing triple-separability in the VPA context |
|figs| Contains general figures and model outputs |
|R_scripts| Contains scripts for model runs, visualizations, and constructing a precision matrix with marginal variance |
|src| Contains source code for TMB model, which is compiled through R|
|test|R scripts testing the development of the TMB model and checks for correct model dimensions|
|output| Model outputs from 2x2x2 factorial model formulations |

# Specification of the Precision Matrix
| Purpose  | Location |
| --------| --------|
|Precision Matrix| Construction of the precision matrix with marginal variance can be [found here](https://github.com/chengmatt/GMRF_WAA/blob/master/R_scripts/make_precision/Construct_precision.R) as a R function, or [here](https://github.com/chengmatt/GMRF_WAA/blob/26b38d6e5dd64bd79d16f1c1d4f3671364bfb93d/src/GMRF_WAA.cpp#L14-L132) as a local function defined in TMB|
|Conditional Variance &Omega; matrix | Code for implementing the &Omega; matrix with condtional variance can be found [here](https://github.com/chengmatt/GMRF_WAA/blob/26b38d6e5dd64bd79d16f1c1d4f3671364bfb93d/src/GMRF_WAA.cpp#L78-L87)|
|Marginal Variance &Omega; matrix | Code for implementing the &Omega; matrix with marginal variance in TMB can be found [here](https://github.com/chengmatt/GMRF_WAA/blob/6ebbeecfe4dd002de66c655d7d687f688cdd1954/src/triple_sep_waa.cpp#L83-L116?plain=1)|

# Examples
## Gulf of Alaska Walleye Pollock (Selectivity)
* Monnahan, C.C., Adams, G.D., Ferriss, B.E., Shotwell, S.K., McKelvey, D.R., McGowan, D.W., 2023. 1. Assessment of the Walleye Pollock Stock in the Gulf of Alaska 128. Accessed from https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2023/GOApollock.pdf

## US Pacific Sardine (Weight at age)
* Kuriyama, P.T., Allen Akselrud, C., Zwolinski, J.P., Hill, K.T., 2024. Assessment of the Pacific sardine resource (Sardinops sagax) in 2024 for U.S. management in 2024-2025. https://doi.org/10.25923/JYW3-YS65

  


