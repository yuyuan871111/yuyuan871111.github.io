---
title: "FastEval Parkinsonism"
date: 2024-09-17T15:18:11Z
draft: False
# github_link: "https://github.com/gurusabarish/hugo-profile"
author: "Yu-Yuan (Stuart) Yang"
tags:
  - Parkinson's disease
  - Artificial intelligence
  - motor assessment
  - FASTEVALP
image: /images/blog_pics/fastevalp.png
description: "Fast evaluation of motor assessment at home."
toc: 
---

The Motor Disorder Society’s Unified Parkinson’s Disease Rating Scale (MDS-UPDRS) is commonly used to assess bradykinesia, a hallmark symptom of Parkinson’s disease (PD). However, it falls short in capturing the full variability of bradykinesia throughout the day outside of a clinical setting. To address this limitation, we present FastEval Parkinsonism ([https://fastevalp.cmdm.tw](https://fastevalp.cmdm.tw)), a deep learning-based video system that enables users to capture key motion points, estimate severity, and generate summary reports. Using 840 finger-tapping videos from 186 participants (103 with Parkinson’s disease (PD), 24 with atypical parkinsonism (APD), 12 elderly individuals with mild parkinsonism signs (MPS), and 47 healthy controls (HCs)), we apply a dilated convolution neural network combined with two data augmentation techniques. Our model achieves acceptable accuracy (AAC) of 88.0% and 81.5%. The frequency-intensity (FI) measure of thumb-index finger distance emerged as a critical hand parameter for quantifying performance. Additionally, our model demonstrated versatility with multi-angle videos, validated using an external database of over 300 PD patients.

## Know more from the reference...
1. **Yang, Y. Y.**, Ho, M. Y., Tai, C. H., Wu, R. M., Kuo, M. C., and Tseng, Y. J. (Feb. 2024). "FastEval Parkinsonism: an instant deep learning–assisted video-based online system for Parkinsonian motor symptom evaluation." NPJ Digital Medicine, DOI: [10.1038/s41746-024-01022-x](https://doi.org/10.1038/s41746-024-01022-x)  
2. **Yang, Y. Y.** (April 2023). “FastEval Parkinsonism: an instant deep learning-based online self-evaluation system for the diagnosis of Parkinson's symptoms with hand videos using finger tapping.” NTU BEBI (***2023 Best Master Thesis Award***), DOI:[10.6342/NTU202300716](https://doi.org/10.6342/NTU202300716)  