---
title: "What is the implicit hydrogen? Why need an additional interaction class? [GSoC 2025 - Week 8]"
date: 2025-07-24T15:18:11Z
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
image: /images/blog_pics/gsoc-wk-8.png
description: "The progress of my GSoC project at week 8"
toc: 
---
Hi again, I am now in the middle of my GSoC journey. Now, we are going to have the main course: **implicit hydrogen bond interactions**.

In this blog post, you will know more about why calculating the implicit hydrogen bond interaction is important and how we design the algorithm of implicit hydrogen bond interactions.


## What is the implicit hydrogen?
From [OpenEye's documentation](https://docs.eyesopen.com/applications/picto/picto/implexplHydro.html), implicit hydrogens are **"the hydrogens that are associated with the parent atom, they are not stored as an atom"**. This is quite common when we store a molecule in a file. We only use heavy atoms to indicate the topology of the molecule. This is a more efficient way to save information (because we can determine the number of hydrogens and their bonding to heavy atoms solely based on the valence and formal charge of the heavy atoms). 

## Why hydrogen bond interation is important?
Hydrogen bond interactions are a specific type of molecular interaction between hydrogen bond donors (D-H) and acceptors (A). This type of molecular interaction is a commonly-used and crucial feature to investigate how the molecule interats with the other moleules. It is also used as a property for a molecule in small-molecule drug design. 

For example, Lipinski's Rule of Five (a famous small-molecule drug discovery guideline) suggests [[ref]](https://www.sciencedirect.com/science/article/pii/S0169409X12002797?casa_token=IdgIgq5ZrIYAAAAA:relfN93f6TpFQfuzx8AQaudyn7mL4GYkE5DpBG5f9YetjEnSJX8a6g1pcXd7OUuaCrAWeIa7G5Y):
> Poor absorption or permeation is more likely when     
  (1) there are more than 5 **H-bond donors**    
  (2) there are more than 10 **H-bond acceptors**    
  (3) the molecular weight (MWT) is greater than 500   
  (4) the calculated Log P (CLogP) is greater than 5 (or MlogP>4.15). 

### Why "implicit" hydrogen bond interaction?
In experimental structure, you don't know where is the actual positions of your hydrogens most of the time. It is because hydrogen's low scattering power (if using X-ray crystallography) and the resolution of the structure (learn more from [here](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/crystallographic-data)). A normal bond length is 1-2 Å and the bond length involving hydrogens is even shorter. If the resolution is not good enough (<2.5 Å), it is hard to determine the actual position of atoms.

Thus, if we need to calculate the hydrogen bond interaction (explicit method) for an experimental structure, we will need to use protonation software (e.g., [PROPKA](https://github.com/jensengroup/propka), [PypKA](https://github.com/mms-fcul/PypKa), etc.) to model the hydrogens.

> To make users' lifes more easier, we think it is important to have "implicit" hydrogen bond interaction in ProLIF. 

Users will be able to compare the experimental structure with computational structure (e.g., from AlphaFold, RoseTTAFold, etc.) without additional data preparation.

## How to get hydrogen bond interactions?
A common method to detect this interaction is to use **atom type** and **geometrical constraints** to detect it. 

For instance, in [ProLIF's implementation (explicit version)](https://prolif.readthedocs.io/en/stable/_modules/prolif/interactions/interactions.html#HBAcceptor), the acceptor and donor are detected first. Only a limited elements can be the acceptor, such as N, O, S, etc (look at the SMARTS pattern for the definition of the acceptor, [visualization](http://smarts.plus/view/cf9f0905-def4-441e-92ec-db1a9af13afa) in below figure). A donor should be O, S or N with at least one hydrogen. 


Then, the geometry between the donor and acceptor pair is checked. A common threshold (which may vary across different software) is that the distance between the donor and acceptor should **no exceed 3.5 Å** and the D-H---A angle should range **between 130° and 180°**. 


### The SMARTS pattern for the acceptor in ProLIF (2.0.0)
![acceptor's SMARTS pattern: [#7&!$([nX3])&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4&+1]),O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O),o+0,F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]](/images/blog_pics/gsoc-wk-8-001.png)

> However, there will be an issue for a molecule without explicit hydrogens.

1. How do you detect your hydrogen bond donor?
2. How to check the geometrical constraints?

**Once these issues fixed, we can detect the hydrogen bond interactions with implicit hydrogens.**


## How do we design implicit hydrogen bond calculations?
Following the same concept of the explicit hydrogen bond calculation, we will need to **detect the donor/acceptor first** and then **check the geometry of the donor-acceptor pair**.

There is no need to modify the acceptor's SMARTS pattern, but the donor's SMARTS pattern must be converted into an implicit version [[visualization]](https://smarts.plus/view/c717e25a-07ef-4f15-9963-02d0ae92a4e5) (also in below figure).

#### New design for H-bond donor
![implicit-H donor's SMARTS pattern: [$([O,S;+0&h]),$([N;v2&h,v3&h,v4&+1&h]),n+0&h]](/images/blog_pics/gsoc-wk-8-002.png)

Now, the first issue solved. Let's deal with the second one.


### Geometrical constraints for implicit hydrogen bonds
The geometrical constraints for "explicit" hydrogen bonds have two parts: **(1) D-A distance**, **(2) D-H---A angle.** Let's see whether we can use the same constraints in the implicit hydrogen bonds calculation. 

#### D-A distance
The calculaiton of D-A distance doesn't involve the hydrogens, so we don't need to modify anything. 

#### D-H---A angle or related geomtrical constraints
However, D-H---A angles are tricky because we don't have a determinstic coordinates of hydrogens now. 

We do know some constraints about the hydrogens. Let's discuss some situations:

##### sp-hybridized donor (linear)
For an sp-hybridized donor, the hydrogen is always aligned along the axis between the donor and its neighboring atom. The X-D-H is 180° (where X is the neighboring atom of the donor). To form a perfect hydrogen bond interaction, the acceptor's lone pair electrons should align with this axis. If the acceptor's lone pair electrons are far way from this axis, the hydrogen bond is impossible to form (that is, ideal X-D-A is 180°).


##### sp<sup>2</sup>-hybridized donor (plane)
For an sp<sup>2</sup>-hybridized donor, its neighboring hydrogens or heavy atoms are always located on the same plane as the donor. The X-D-H is 120° (where X is the neighboring atom of the donor). Thus, the ideal hydrogen bond interaction, the acceptor's lone pair electrons should be on the plane and along with the D-H axis (that is, ideal X-D-A is 120°). 


##### sp<sup>3</sup>-hybridized donor (tetrahedral)
For an sp<sup>3</sup>-hybridized donor, its neighboring hydrogens or heavy atoms are always located on the same plane as the donor. The X-D-H is 109.5° (where X is the neighboring atom of the donor). Thus, the ideal hydrogen bond interaction, the acceptor's lone pair electrons should be align with the D-H axis (that is, ideal X-D---A is 109.5°). However, if there are two or more hydrogens on the donor, multiple D-H axses exist (but the ideal X-D-A maintains 109.5°).

#### X-D---A angle helps geometrical checks (replacing D-H---A angle)
As discussed in the previous sessions, different hybridized donor could have different ideal X-D-A angles. If we compare the ideal X-D---A angle and the actual one, we will know whether it is possible to form hydrogen bond interaction (from the pointview of the donor atom). 
![sp/sp2/sp3-hybridized donor](/images/blog_pics/gsoc-wk-8-003.jpg)

> In a nutshell, setting a tolerance for deviations between ideal and actual geometric properties helps exclude unlikely hydrogen bond conditions.

#### To be more generalized and in a more mathmatical way:
We calculated **"atom angles"** and **"plane angles"** (called as "geometrical properties") of the donor and acceptor. The concept of "atom angles" and "plane angles" are come from [Mol*](https://molstar.org/), and they represent the actual geometrical relationship between the donor and acceptor. The deviation between the ideal and actual "geometrical properties", determined by taking the absolute difference, can be an indicator to evaluate the possibilty for forming a hydrogen bond (without explicit hydrogens).

**The physical meaning of the deviation of each geometrical properties:**

* **Donor atom angles (X-D---A angle)** -> evaluate how closely the actual **hydrogen donation** (partial positive) aligns with the ideal geometry (ideal sp-hybridized donor atom: 180°, sp<sup>2</sup>: 120°, sp<sup>3</sup>: 109.5°).
* **Acceptor atom angles (D---A-X angle)** -> evaluate how closely the actual **lone pair electrons acceptance** (partial negative) aligns with the ideal geometry (ideal sp-hybridized acceptor atom: 180°, sp<sup>2</sup>: 120°, sp<sup>3</sup>: 109.5°).
* **Donor plane angles (sp<sup>2</sup>-hybridized donor only)** -> evaluate whether the **hydrogen of the donor** lies on the same plane as the acceptor.
* **Acceptor plane angles (sp<sup>2</sup>-hybridized acceptor only)** -> evaluate whether the **lone pair electrons of the acceptor** lies on the same plane as the donor.

#### Donor atom angle and its deviation (2D diagram) 
![donor atom angle](/images/blog_pics/gsoc-wk-8-004.png)

where X is the neighboring atom of the donor atom.

#### Donor plane angle (3D diagram) 
![donor plane angle](/images/blog_pics/gsoc-wk-8-005.png)

where X1 and X2 are the neighboring atoms of the donor atom.

## Comparison between the explicit and implicit H-bond interaction classes in ProLIF
Here is an example of H-bond interaction detection between the explicit and implicit methods. You will notice that no hydrogens are displayed for the ligand and proteins on the right-hand side (they are implicit!!!). **Also, you can clearly see there are more interactions found by the implicit H-bond method.**

However, if you look into the H-bond donors of the ligand forming additional implicit H-bond interactions (the two oxygens on the right), you will find they are both sp<sup>3</sup>-hybridized donors, which suggests the hydrogen of the donors can be rotatable during molecular dynamics simulation (or said, having multiple different conformational states of the hydrogen). **Considering an alternative conformational state of the hydrogen, it is also possible to form a hydrogen bond interaction, as suggested by the implicit method.** Under these circumstances, the implicit H-bond method could somehow provide more information (by mimicking multiple states) than the explicit H-bond method (by focusing on single state).

![Implicit H-bond interaction in ProLIF](/images/blog_pics/gsoc-wk-8-006.png)

Note: We leave an option for users to include or exclude water molecules. This example includes water molecules for the implicit method but excludes them for the explicit method.

### More than donor/acceptor detection and geometrical checks for the implicit H-bond method
To reduce the false positive of the calculation, we include a feature (inspired by Autodock Vina) for the detected implicit H-bond interaction.

In [Autodock Vina](https://vina.scripps.edu/), H-bond interactions are considered when calculating the docking score. Their H-bond term uses the D-A distance and the Van der Waals radii of donors and acceptors to interpolate the hydrogen bond contribution to the docking score. This concept can also serve as an indicator to assess the "potential" of H-bond interactions in our implicit H-bond method.

![Autodock Vina: H-bond terms](/images/blog_pics/gsoc-wk-8-007.png)


## Summary
Now, we have completed the core development of implicit H-bond interaction method. We will next focus on validating this method with explicit hydrogen bonds calculation (and tuning some parameters in the implicit method, such as tolerance, cutoff of the Autodock Vina H-bond potential term).


***See you around.***

## More about GSoC 2025: MDAnalysis x ProLIF project 5
See how we manage this project: click on this [[link]](https://github.com/yuyuan871111/GSoC2025_Hbond_PM).

Know more about the project details via [GSoC Project Page](https://summerofcode.withgoogle.com/programs/2025/projects/5Otkx8vp) and [ProLIF Github](https://github.com/chemosim-lab/ProLIF).


Note: The top picture was generated by **DALL.E3** with **Watercolor** style, aspect ratio: **16:9**, generation count **4** (I picked one picutre out of four.), and prompts: "I would like a drawing with the idea of the below keywords but no texts: implicit hydrogen versus explicit hydrogen, geometrical constraints, coding for cheminformatics, molecule". After the first generation, the picture was refined by removing manually-brushed area and upscaling.


## Last and next posts of GSoC
* Previous: [week 6](/blogs/gsoc-week-6)
* Next: [week 10](/blogs/gsoc-week-10)