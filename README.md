data and code of Pfitznerl et. al. 2026

- code repository, for the first part without the SEF
  
    the code largely uses the opera package code of Opera version 1.1.1 (Gaillard et. al. 2016).
    
    The BOA function is different from the BOA function of opera.
    We have implemented the SOCO BOA of Wintenberger 2024 whereas in opera, it is the BOA of Wintenberger 2017 which is implemented.
    
    agreg_fenetre_glissante.R is a main script that can run several expert aggregations over several stations and lead times, it also computes different scores.
    
    fonction_prevision_MODELE.R runs the MODELE.R functions and enables to deal with sliding windows, where MODELE can be BOA, EWA, MLprod...
    plot_agreg.R contains plot functions used for the paper
    
    Opera_dev.R is a copy of the opera code with some changes.


- code_SEF repository, for the second part with the SEF
  
    the code largely uses the opera package code of Opera version 1.1.1 (Gaillard et. al. 2016).
    
    The BOA function is different from the BOA function of opera.
    We have implemented the SOCO BOA of Wintenberger 2024 whereas in opera, it is the BOA of Wintenberger 2017 which is implemented.
    We also adapted the function to wake up sleeping experts with transition_leaning_xgboost and probtransition objects.
    
    agreg_fenetre_glissante.R is a main script that can run several expert aggregations over several stations and lead times, it also computes different scores.

    transition_learning_xgboost_multiple_ech_sta to train xgboost with already available training data.
    transition_learning_xgboost_multiple to train xgboost with data gatered online, only from the same staion and lead time.

    To train xgboost only with data from the same station and lead times comment £ and us the first BOA function.
    To train xgboost with more data, that is already in the right format, uncomment £ and use the second BOA funciton.
    
    fonction_prevision_BOA.R runs the BOA function and enables to deal with sliding windows.
    plot_agreg.R contains plot functions used for the paper

    test.R, for Diebold Mariano tests and quantile tests
    
    Opera_dev.R is a copy of the opera code with some changes.


