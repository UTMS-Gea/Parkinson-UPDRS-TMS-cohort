# Parkinson's Disease Cohort Study: Clinical and TMS Data Analysis

This repository contains a database with clinical and transcranial magnetic stimulation (TMS) measurements, along with the R code used for statistical analysis of a cohort of patients with Parkinson's disease (PD). The data and analyses are associated with the article:

> **"Longitudinal Dynamics of Clinical and Neurophysiological Changes in Parkinson's Disease: A 4.5-Year Cohort Study"**

## Study Overview
This longitudinal observational study followed a cohort of PD patients:
- **Retrospectively** from 2018 to 2021
- **Prospectively** from 2021 to 2023

The primary objective was to evaluate clinical and neurophysiological parameters over time, with key assessments including:
- **MDS-UPDRS** (Movement Disorder Society-Unified Parkinson’s Disease Rating Scale)
- **TMS-derived neurophysiological measures**

The study consisted of four distinct assessment periods:
- **T1**: September 2018 – June 2019
- **T2**: June 2019 – March 2020
- **T4**: March 2021 – March 2022
- **T5**: March 2022 – March 2023

> **Note**: In-person assessments were suspended from March 2020 to March 2021 due to COVID-19 restrictions.

For more details, please refer to the published article.

---
## Database Structure
The dataset consists of three CSV files located in the `data_UPDRS_TMS` folder.

### 1. `UPDRS.csv`
Contains MDS-UPDRS measurements with 13 columns:
| Column | Description |
|--------|-------------|
| `ID` | Subject ID (22 subjects) |
| `Period` | Evaluation period (T1, T2, T4, T5) |
| `PeriodNum` | Coded numerical representation of the evaluation period |
| `Sex` | Biological sex (M: male, F: female) |
| `MAS` | More affected side |
| `HY` | Hoehn and Yahr scale |
| `P1` | MDS-UPDRS Part I score |
| `P2` | MDS-UPDRS Part II score |
| `P3` | MDS-UPDRS Part III score |
| `P3UL` | Part III score in upper limbs |
| `P3LL` | Part III score in lower limbs |
| `P3right` | Part III score on the right side |
| `P3left` | Part III score on the left side |

### 2. `TMS.csv`
Contains TMS measurements with 15 columns. The first five columns are identical to `UPDRS.csv`, while the remaining columns include:

| Column | Description |
|--------|-------------|
| `mtMA` | Resting motor threshold (more affected hemisphere, % of max stimulator output) |
| `mtLA` | Resting motor threshold (less affected hemisphere, % of max stimulator output) |
| `latMA` | Latency of MEP (more affected hemisphere, ms) |
| `latLA` | Latency of MEP (less affected hemisphere, ms) |
| `durMA` | Duration of MEP (more affected hemisphere, ms) |
| `durLA` | Duration of MEP (less affected hemisphere, ms) |
| `silMA` | Cortical silent period (more affected hemisphere, ms) |
| `silLA` | Cortical silent period (less affected hemisphere, ms) |
| `aucMA` | Area under the recruitment curve (more affected hemisphere, mV.ms) |
| `aucLA` | Area under the recruitment curve (less affected hemisphere, mV.ms) |

> **MEP**: Motor Evoked Potential

### 3. `IOCurves.csv`
Contains MEP amplitude data from recruitment curves (input-output curves) with 5 columns:

| Column | Description |
|--------|-------------|
| `ID` | Subject ID |
| `Period` | Evaluation period (T1, T2, T4, T5) |
| `Intensity` | Stimulation intensity (% of resting motor threshold) |
| `ampMA` | Amplitude in the more affected hemisphere (mV) |
| `ampLA` | Amplitude in the less affected hemisphere (mV) |

---
## R Code for Data Analysis
The R code used for analysis consists of four scripts:

| File | Description |
|------|-------------|
| `EP_cohort_UPDRS.r` | Analyzes MDS-UPDRS data and generates Figure 1 and Supplementary Figure 1 |
| `EP_cohort_TMS.r` | Analyzes TMS data and generates Figure 2 and Supplementary Figures 2, 3, and 4 |
| `EP_cohort_IOcurves.r` | Generates Figure 3 related to recruitment curves |
| `EP_cohort_correlations.r` | Analyzes within- and between-subject correlations between clinical and TMS variables |

### **Execution Order**
To ensure proper execution, run the scripts in the following order:

1. `EP_cohort_UPDRS.r`
2. `EP_cohort_IOcurves.r`
3. `EP_cohort_TMS.r`
4. `EP_cohort_correlations.r`

This sequence ensures that all required libraries are loaded and that necessary data frames created in one script are available for subsequent analyses.

---
## License
MIT License (For R Scripts)
CC-BY-NC 4.0 License (For Dataset)

## Citation
If you use this dataset or code, please cite the corresponding article:
> [Citation pending.]

---
## Contact
For questions or clarifications, please contact:

> *Oscar Arias-Carrion*
ariasemc2@gmail.com
Hospital General Dr. Manuel Gea González

> *Emmanuel Ortega-Robles*
emmanuel.ortega@salud.gob.mx
Hospital General Dr. Manuel Gea González

---
### Acknowledgments
We acknowledge the contributions of all study participants and research team members involved in data collection and analysis.

