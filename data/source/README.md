

README
======


This file describes the content of the different files included in this repository to 
reproduce results from [1] and some of its supplementary results.

## Input files ##

* 20151016\_Functions\_remainder.csv

	Functions measured in [2]. The relevant quantities used in [1] are labelled with "7", and include:
	
	* Community: Id of the sample
	* Replicate
	* Plate
	* mgCO2.7: CO2 measured along 7 days of experiment
	* CPM7: Cell counts at the end of the experiment 
	* pgRPC.7: CO2 per cell
	* ATP7: ATP measured (nM)
	* mG7: beta glucosidase (mM)
	* mN7: beta chitinase (mM)
	* mX7: xylosidase (mM)
	* mP7: phosphatase (mM)
  
* samples\_metadata\_time0.tsv

	* Samples: Id of the sample	
	* Part.dates: Date of sampling
	* Part.GPS.PAM: Optimal sampling sites	
	* Part.SparCC.PAM.t0: Optimal partition using SparCC 
	* Part.SJD.PAM.t0: Optimal partition using Jensen-Shannon Divergence
	* Part.Dir.t0: Optimal partition using Dirichlet mixtures
	* Part.month: Month in which the community was sampled
* Dist\_GPS-Haversine.dat

	Haversine (spatial) distances between samples

* corMat-SparCC\_20151016\_OTU\_remainder.clean.samples.txt

	Matrix of correlations between samples computed with SparCC
	
* distMat\_ShannonJensen\_Samples\_Time0.clean.dat

	Distance matrix computed with Jensen-Shannon divergence.


## Supplementary results ##

* SEMmodels.zip

	Results for the Structural Equation Models analysed. The structure of the folders follows the one 
	available at the repository of the [project ](https://github.com/apascualgarcia/TreeHoles_descriptive).
	
* TaxaSummaries.zip

  The file contains one folder for each community-class, with matrices in different formats (biom and txt)
  computing the relative abundances of the OTUs at different taxonomic levels (labelled L2 being the proxy for
  Phylum to L6, the proxy of species). These matrices can be visualized interactively opening with a web
  browser the files area_charts.html.

#### References ####

> [1] Pascual-García, A., & Bell, T. (2019). Community-level signatures of ecological succession 
>    in natural bacterial communities. bioRxiv, 636233

> [2] Rivett, Damian W., and Thomas Bell. Abundance determines the functional role of bacterial 
>     phylotypes in complex communities." Nature microbiology 3.7 (2018): 767.
