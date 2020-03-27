# Bayesian Factor Analysis for Inference on Interactions
## Federico Ferrari, David B Dunson

This article is motivated by the problem of inference on interactions among chemical exposures impacting human health outcomes. Chemicals often co-occur in the environment or in synthetic mixtures and as a result exposure levels can be highly correlated. We propose a latent factor joint model, which includes shared factors in both the predictor and response components while assuming conditional independence. By including a quadratic regression in the latent variables in the response component, we induce flexible dimension reduction in characterizing main effects and interactions. We propose a Bayesian approach to inference under this Factor analysis for INteractions (FIN) framework. Through appropriate modifications of the factor modeling structure, FIN can accommodate higher order interactions and multivariate outcomes. We provide theory on posterior consistency and the impact of misspecifying the number of factors. We evaluate the performance using a simulation study and data from the National Health and Nutrition Examination Survey (NHANES). Code is available on GitHub.

[Link to manuscript](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2020.1745813#.Xn4EjS2ZNQI)

*A faster implementation of the algorithm is available at https://github.com/poworoznek/infinitefactor and will be soon be available on CRAN*
