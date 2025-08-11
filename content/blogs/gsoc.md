---
title: "MDAnalysis x ProLIF Project 5: H-Bond interactions from implicit hydrogens"
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
**"MDAnalysis x ProLIF Project 5: H-Bond interactions from implicit hydrogens"**. It is a cross-link project between MDAnalysis (an open source organization, also a python package) and ProLIF (a python package).

### What is ProLIF? What I am doing in my GSoC project?
ProLIF, a tool for analyzing protein-molecule interactions in molecular dynamics, struggles with detecting hydrogen bonds when only heavy atom data is available. Many users compare structures without explicit hydrogens, but current protonation tools are complex and error-prone for beginners. This proposal introduces an **"implicit hydrogen bond detection"** method using heavy atom positions, balancing efficiency and accuracy. Inspired by tools like [Mol*](https://molstar.org/) and [Chimera](https://www.cgl.ucsf.edu/chimera/), it will use distance-based calculations with geometrical constraints and a confidence score inspired by [AutoDock Vina](https://vina.scripps.edu/). To improve usability, a "Protein Helper" function will standardize residue names and assign protonation states from common force fields. The method will be validated using PINDER and PLINDER datasets against standard H-bond detection methods.


### Know more about our GSoC project...
1. **See how we manage this project:** click [here](https://github.com/yuyuan871111/GSoC2025_Hbond_PM).
2. [GSoC Project Page](https://summerofcode.withgoogle.com/programs/2025/projects/5Otkx8vp): to know the proposal and the project details.
3. [ProLIF Github](https://github.com/chemosim-lab/ProLIF): main codebase. Feel free to raise *Issue* (it would help ProLIF making better.) and contribute codes with *Pull Request*.


### GSoC blog posts (up to date)
1. [Week 0~4]: [First few weeks](/blogs/gsoc-week-4) of my GSoC project.   
2. [Week 5~6]: What happen at/until [week 6](/blogs/gsoc-week-6)?
3. [Week 7~8]: Next stage of the project: build more at/until [week 8](/blogs/gsoc-week-8).
4. [Week 9~10]: Complete and validate my method at/until [week 10](/blogs/gsoc-week-10).
5. [Week 11~12]: TBU [week 12](/blogs/gsoc-week-12)
6. [At the end of GSoC]: TBU [GSoC summary](/blogs/gsoc-week-summary)