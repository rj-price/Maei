# InterProScan Annotation for Predicted Effectors

InterProScan Analysis on the set of 90 non-redundant predicted effectors from Maei and FoEC. 

The combined set of predicted effectors from Maei (13.08.2020) and FoEC (04.08.2020) (see -i flag in command) were searched using InterProScan default paramenters. 

The -goterms flag was used for GO analysis to try and predict function. 

Comman (on Vettel):

[interproscan.sh](http://interproscan.sh/) -i /home/u1983390/Fusarium_data/MIMPS/MIMP_Maei_Searches_13.08.2020/Compare_Maei_to_FoEC/Original_Clusters/all_Maei_and_FoEC_comb.NON-REDUNDANT.CLUSTERED.fastaÂ  -goterms -b ./

---

**OutputÂ** 

Output can be found in the following directory;

Vettel:

/home/u1983390/Fusarium_data/MIMPS/FoEC_and_Maei_downstream/InterProScan
