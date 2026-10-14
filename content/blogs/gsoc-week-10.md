---
title: "Validation and Optimization: tuning parameters for the implicit H-bond interaction method using datasets [GSoC 2025 - Week 10]"
date: 2025-08-08T08:08:08Z
draft: False
github_link: "https://github.com/chemosim-lab/ProLIF"
author: "Yu-Yuan (Stuart) Yang"
tags:
  - Google Summer of Code 2025
  - MDAnalysis
  - ProLIF
  - H-bond interaction
  - GSOCMDAPL
  - implicit hydrogen
  - parameter tuning
image: /images/blog_pics/gsoc-wk-10.png
description: "The progress of my GSoC project at week 10"
toc: 
---
Hello, welcome to blog post on the week 10 of my GSoC journey! From the [last blog post](/blogs/gsoc-week-8), we completed the development of the implicit H-bond interaction class. However, it is not perfect. Before getting this method merged into the main ProLIF, we have to validate our method using the explicit H-bond method. Thus we used a dataset, [PLINDER](https://www.plinder.sh/), to achieve this.

In this blog post, you will learn how we validated our implicit H-bond method and understand what parameters were tuned.

## How to validate our method?
We used explicit H-bond interaction method as ground truth to validate our method using a protein-ligand complex. Here are the steps:
1. Protonate the protein and add hydrogens for the ligand.  
2. Calculate explicit H-bond interactions.
3. Remove the hydrogens in the ligand and protein, and then calculate implicit H-bond interactions (our method).
4. Compare donor-acceptor pairs from explicit and implicit H-bond methods.

### An example
First, for protein protonation, we used [PDB2PQR](https://pdb2pqr.readthedocs.io/en/latest/) (PROPKA3) and [reduce](https://github.com/rlabduke/reduce) (if PDB2PQR is not working for the system). For adding hydrogens for ligands, I used [openbabel](https://openbabel.org/). See how we dealt with a lot of data: [click here](https://github.com/yuyuan871111/GSoC2025_Hbond_PM/blob/main/validation/utils/protonate.py). 

Then, we read the protein with `RDKit`, converted it into `prolif.Molecule`, and fixed bond orders with our `ProteinHelper`. The ligand was also read as a `prolif.Molecule`.

![read protein](/images/blog_pics/gsoc-wk-10-001.png)

![read ligand](/images/blog_pics/gsoc-wk-10-002.png)
It is noted that multiple ligands can interat with a single protein. We used a list to save all possible ligands for that protein receptor.

Next, we calculated the explicit H-bond interations with default parameters (distance = 3.5 Å, D-H-A angle is between 130° and 180°). They were considered as ground-truth H-bond interactions for this protein-ligand system.  
![explicit H-bond](/images/blog_pics/gsoc-wk-10-003.png)

The above results show the interaction pairs between the ligand and residues. 

Then, we can also do the same thing for implicit H-bond method. We used RDKit to remove the explicit hydrogens before `ProteinHelper`.
![read protein for implicit Hbond](/images/blog_pics/gsoc-wk-10-004.png)

Same for the ligand, we removed the explicit hydrogens.
![read ligand for implicit Hbond](/images/blog_pics/gsoc-wk-10-005.png)

For the calculation of the implicit H-bond interations, we included water molecules. Also, we set a relatively large tolerance for each deviation of atom/plane angles so that we wouldn't miss any possible pairs from the explicit methods (yes, we intentionally minimized false negatives but it might generate lots of false positives).

![implicit H-bond](/images/blog_pics/gsoc-wk-10-006.png)

Here is the list of interaction pairs from implicit H-bond method.

![implicit H-bond2](/images/blog_pics/gsoc-wk-10-007.png)

Finally, we compared the different sets of interaction pairs from implicit and explicit methods. In implicit H-bond method, 23 interactions are detected. However, only 13 interactions are actually detected by the explicit H-bond method.

![confusion matrix](/images/blog_pics/gsoc-wk-10-008.png)

We also calculated **"Tanimoto coefficient"** (or said Jaccard index, or simply Intersection over Union (IoU)) between two methods.

```python
Tanimoto = 13 / (13 + 10 + 0) ~= 0.565
``` 

**Now, we had a quatitative method to compare interaction pairs from two differnt methods** (not only maunally visulization like the below picture).

![interaction by visual](/images/blog_pics/gsoc-wk-10-009.png)


## Tuning the implicit H-bond parameters
The next thing to think about is: 

> How do we improve the accuracy of our method?

In the implicit H-bond methods, we have several parameters, including **acceptor atom angle deviation, donor atom angle deviation, acceptor plane angle, donor plane angle**. If we know how to set a good threshold for each parameter, we can have more accurate results. To achieve this, we would like use the power of datasets.

**Our goal is to analyze all hydrogen bond interaction pairs** identified using the explicit hydrogen bond method in the dataset (we used [PLINDER](https://www.plinder.sh/) dataset). Specifically, we examined the distribution of geometric parameters (**acceptor atom angle deviation, donor atom angle deviation, acceptor plane angle, donor plane angle**). Additionally, we compared these "ground-truth" H-bond pairs with false positives generated by the implicit H-bond method to determine whether their distribution profiles differ significantly. In the end, we can decide cutoff values for each geometric parameters (that is why we called parameter tuning).


### PLINDER dataset
[PLINDER](https://www.plinder.sh/) is a protein-ligand dataset with mutliple annotations. The dataset was splitted into training (n~=420K), validation (n=1157) and testing (n=1436) sets to minimize leakage and maximize the quality of test. 

For our parameter tuning, we used PLINDER's validation set as our first dataset to investigate the parameters' distribution, and PLINDER's testing set for testing the selection of the parameters.

### Donor and acceptor atom angle deviation's distribution
The points represent each interaction pair's acceptor atom angle deviation and donor atom angle deviation in the PLINDER's validation set (overall averaged Tanimoto coefficient: 0.565). The blue one is detected by both implicit and explicit H-bond methods while the red one is detected by the implicit H-bond method only. From the distribution of the acceptor atom angle deviation, two sets are identical. This suggests that we cannot select a cutoff of this parameter to reduce the false positives. However, we can select **a cutoff value of 25°** for the donor atom angle deviation. Because lots of false positives have a greater donor atom angle deviation compared to the ground truth.

![atom angle deviation](/images/blog_pics/gsoc-wk-10-010.png)
It is noted that the tolerance of atom angle deviation was originally set in 60°. There are 14 (~1%) protein-ligand system fails to detect H-bond interactions wither from the explicit or implicit methods due to computing errors (under investigation).

### Donor and acceptor plane angles' distribution
We did a similar analyses for plane angles. Interesting, **the acceptor plane angles again show no difference between the ground truths and the false positives.** But, we can set **a donor plane angle of 30°** to reduce false positives.

![plane angles](/images/blog_pics/gsoc-wk-10-011.png)
It is noted that the tolerance of acceptor and donor plane angles was originally set in 90° and 60°, respectively. There are 14 (~1%) protein-ligand system fails to detect H-bond interactions wither from the explicit or implicit methods due to computing errors (under investigation).

## Summary
Now, we have found the cutoff values for some parameters in the implicit H-bond interaction method. 
* Tolerance of acceptor atom angle deviation: it doesn't matter.
* Tolerance of donor atom angle deviation: 25°.
* Tolerance of acceptor plane angle: it doesn't matter.
* Tolerance of donor plane angle: 30°.

We will next finalize our functions and it will then be ready to go.


***See you in the next blog post.***

## More about GSoC 2025: MDAnalysis x ProLIF project 5
See how we manage this project: click on this [[link]](https://github.com/yuyuan871111/GSoC2025_Hbond_PM). The code for the validation shown in the blog post is also in the same repo. If you want to know more about the project's details, you can find it via [GSoC Project Page](https://summerofcode.withgoogle.com/programs/2025/projects/5Otkx8vp) or [ProLIF Github](https://github.com/chemosim-lab/ProLIF/tree/gsoc_implicit_hbond).


Note: The top picture was generated by **Flux Pro 1.1 Ultra** with **Oil Painting** style, aspect ratio: **16:9**, generation count **4** (I picked one picutre out of four.), and prompts: "I would like a drawing with the idea of the below keywords but no texts: implicit hydrogen versus explicit hydrogen, geometrical constraints, coding for cheminformatics, molecule, parameters tuning, using datasets, Protein-Ligand interaction, Protonation tools.".


## Last and next posts of GSoC
* Previous: [week 8](/blogs/gsoc-week-8)
* Next: [week 12](/blogs/gsoc-week-12)