---
title: "[Google Summer of Code 2025] MDAnalysis x ProLIF Project 5: H-Bond interactions from implicit hydrogens"
date: 2025-05-18T15:18:11Z
draft: False
github_link: "https://github.com/chemosim-lab/ProLIF"
author: "Yu-Yuan (Stuart) Yang"
tags:
  - Google Summer of Code 2025
  - MDAnalysis
  - ProLIF
  - H-bond interaction
  - GSOCMDAPL
image: /images/blog_pics/gsoc.png
description: "Google Summer of Code 2025 with MDAnalysis and ProLIF"
toc: 
---

Google Summer of Code is a global program that connects contributors with open source organizations for 12+ week coding projects guided by mentors.

This year (2025), I am glad to be a contributor in 
**"MDAnalysis x ProLIF Project 5: H-Bond interactions from implicit hydrogens"**. 

ProLIF, a tool for analyzing protein-molecule interactions in molecular dynamics, struggles with detecting hydrogen bonds when only heavy atom data is available. Many users compare structures without explicit hydrogens, but current protonation tools are complex and error-prone for beginners. This proposal introduces an **"implicit hydrogen bond detection"** method using heavy atom positions, balancing efficiency and accuracy. Inspired by tools like [Mol*](https://molstar.org/) and [Chimera](https://www.cgl.ucsf.edu/chimera/), it will use distance-based calculations with geometrical constraints and a confidence score inspired by [AutoDock Vina](https://vina.scripps.edu/). To improve usability, a "Protein Helper" function will standardize residue names and assign protonation states from common force fields. The method will be validated using PINDER and PLINDER datasets against standard H-bond detection methods.

See how we manage this project: click on this [[link]](https://github.com/yuyuan871111/GSoC2025_Hbond_PM).


Know more from the reference: 
1. [GSoC Project Page](https://summerofcode.withgoogle.com/programs/2025/projects/5Otkx8vp).
2. [ProLIF Github](https://github.com/chemosim-lab/ProLIF).