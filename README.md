data and code of Pfitznerl et. al. 2026

the code largely uses the opera package code of Opera version 1.1.1 (Gaillard et. al. 2016).

The BOA function is different from the BOA function of opera.
We have implemented the SOCO BOA of Wintenberger 2024 whereas in opera, it is the BOA of Wintenberger 2017 which is implemented.

agreg_fenetre_glissante.R is a main script that can run several expert aggregations over several stations and lead times, it also computes different scores.

fonction_prevision_MODELE.R runs the MODELE.R functions and enables to deal with sliding windows, where MODELE can be BOA, EWA, MLprod...
plot_agreg.R contains plot functions used for the paper

Opera_dev.R is a copy of the opera code with some changes.


