# ClimaticResponses_Adapations_and_RangePosition

Main code used for the MS: "Species responses to weather anomalies depend on local adaptation and range position" 
Yolanda Melero, Luke C. Evans, Mikko Kuussaari, Reto Schmucki, Constantí Stefanescu, David B. Roy, Tom H. Oliver

Note: The datasets used for this study are available from the European Butterfly Monitor Scheme (eBMS) via a signed license agreement (https://butterfly-monitoring.net/). The dataset generated for the analyses of the study is available via Zenodo at 10.5281/zenodo.15012265. Climatic data are available via ECAD website (https://www.ecad.eu/). For the phylogenetic models, the Butterflies_Europe_tree.nwk is available from Dapporto et al (2019) Doi 10.1111/1755-0998.13059, and from Wiemers et al (2020) doi 10.3897/zookeys.938.50878.

Note: The code includes the general template codes to be used for the analyses:
1. Sims_Fis1.R: Code to create Figure 1, including simulation code for H1 (Fig1. panel b)
2. Bioclimatic_Range.R: Code for the bioclimatic range construction. Data available from eBMS license agreement.
3. Models_H1.R: Models for population growth rate in relation to the local climatic anomalies and range position. Models were done per species category based on the provided code. 
4. Models_2.R: Models for abundance over time in relation to the local climatic anomalies and range position. Models were done per species category based on the provided code.
5. Models_H1_Pheno.R: Models for population growth rate in relation to the local climatic anomalies and range position, including phylogeny (as we observed an effect in the species responses in terms of population growth rate to climcatic anomalies at Melero et al (2022) doi: 10.1038/s42003-022-03088-3. Models were done per species category based on the provided code. 
The code to calculate the index of abundances is available at: https://github.com/RetoSchmucki/rbms
