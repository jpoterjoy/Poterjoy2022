This data assimilation package runs the Local PF with iterations and hybridization with an ensemble Kalman filter, as introduced in Poterjoy (2022).

This version can operate with three different models: 

      Lorenz (1963)
      Lorenz (1996)
      Lorenz (2004)

To run: 

      Open and edit options in main_dyn.m

Directory structure:

      FILTERS <--- contains functions for running EnKF and LPF
      MODELS <--- contains functions for running models
      MISC <--- contains various functions needed by main_dyn.m or filters
      DATA <--- output from experiments are stored here
      FIGS <--- figs from experiments are stored here

