# COVID-19: Who Should Get Vaccinated First? 
**B.Sc. Thesis in Mathematical Statistics | Stockholm University (2021)**

## Advisor
* **Dr. Pieter Trapman**, Stockholm University

## Overview
This repository contains the data and R implementation for my Bachelor’s thesis. The project investigates optimal vaccination strategies for the Stockholm region by modeling disease dynamics across different age and sex cohorts.

* **Model:** An augmented **SEIR (Susceptible-Exposed-Infectious-Recovered)** compartmental model where the standard Recovered state is split into three states: Recovered, Long-term Ill, and Deceased.
* **Objective:** Optimization of vaccination scheme based on three criteria:
  1. Minimizing total transmission.
  2. Protecting vulnerable individuals (long-term illness, older).
  3. Reducing fatality rates.
* **Mathematics:** The final epidemic sizes for each demographic group are computed by numerically solving balance equations.
* **Results:** The optimal vaccination scheme depends on the end goal. The vaccination scheme prioritizing younger individuals reduces transmission the most, while the scheme prioritizing vulnerable individuals protects risk groups the most. However, vaccinating individuals with a higher risk of early infection sufficiently reduces both transmission and protects risk groups.

## Contents
* `2021_16_report.pdf`: The full thesis text (English).
* `/code`: R script for the numerical solvers and data visualization.
* `/Data`: Demographic and epidemiological data from the Stockholm region.
