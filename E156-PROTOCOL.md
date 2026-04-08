# E156 Protocol — PredModelMA

**Project**: PredModelMA
**Created**: 2026-04-09
**Type**: methodological
**Estimand**: Pooled c-statistic (logit transform) and O:E ratio (log transform)
**Dashboard**: https://mahmood726-cyber.github.io/PredModelMA/

## E156 Body (153 words, 7 sentences)

Can prediction model validation studies be synthesized in a browser with the same statistical rigor as R packages like metamisc? We implemented c-statistic pooling via logit transformation (logit(c) = log(c/(1-c))) with delta-method standard errors (SE(logit) = SE(c)/(c(1-c))), and O:E ratio pooling on the natural log scale, both using REML or DerSimonian-Laird random effects. PROBAST risk-of-bias assessment covers four domains (Participants, Predictors, Outcome, Analysis) with traffic-light visualization and domain-level bar charts. Applied to 10 Framingham Risk Score external validations, the tool yields pooled c = 0.742 (95% CI 0.731-0.753) and pooled O:E = 1.09 (1.03-1.16), indicating moderate discrimination with slight underprediction — consistent with published systematic reviews. Forest plots display study-level estimates with weight-proportional symbols and pooled diamonds for both discrimination and calibration. The tool enables rapid prediction model evidence synthesis without software installation or programming. This implementation does not support calibration-in-the-large meta-regression or net benefit decision curve pooling.
