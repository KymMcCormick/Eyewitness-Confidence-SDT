# Eyewitness Confidence under Signal Detection Theory

R code and analyses investigating how eyewitness identification confidence behaves under signal detection theory (SDT) models.

This repository examines how memory strength, evidence variability, decision criteria, lineup size and baseline expectations about culprit presence influence identification decisions and the confidence–accuracy relationship. It includes theoretical probability calculations, simulations and empirical model fitting.

## Research questions

- How do memory strength and baseline expectations about culprit presence affect the frequency and accuracy of high-confidence identifications?
- How do fixed and varying confidence criteria affect culprit and innocent-suspect identification probabilities?
- How do equal-variance and unequal-variance SDT models account for observed identification and confidence distributions?
- How do memory strength and baseline expectations, separately and jointly, influence the mapping of evidence onto identification decisions and confidence?

## Studies

### Confidence mixture model critique

The first project examines whether declining high-confidence positive predictive value (PPV) requires an additional confidence-generating process.

Analyses compare a fixed high-confidence criterion with criteria adjusted to maintain 95% PPV. Empirical analyses examine estimated high-confidence criteria across memory and expectation conditions under equal-variance and unequal-variance SDT models.

**Associated manuscript:** *When confident eyewitnesses are mistaken: is a mixture model necessary?*

Further studies will investigate confidence under different model assumptions and parameter settings.

## Interpretation

PPV is the proportion of suspect identifications that identify the culprit. It depends on both identification probabilities and culprit-present versus culprit-absent lineup base rates.

The probability of identifying an innocent suspect is a separate quantity. Analyses distinguish these outcomes and state the assumptions under which predictions are calculated.

## Reproducibility

Study-specific instructions will describe required R packages, data inputs and how to reproduce analyses, tables and figures. Empirical data availability and access conditions will be documented separately.

## Status

This repository is under development. Code, documentation and analyses may change as projects progress.

## Citation

If you use this code, please cite the repository and the associated publication where applicable. Publication details and version-specific citation information will be added as they become available.

## Author

Dr Kym Michelle McCormick

Adelaide University, South Australia

## Licence

Code is released under the MIT License. See [LICENSE](LICENSE) for the full terms. Third-party materials and any datasets retain their own applicable terms.
