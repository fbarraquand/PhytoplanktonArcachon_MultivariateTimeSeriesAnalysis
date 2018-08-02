### MARSS results

This repository stores the MARSS objects corresponding to our two study sites (Teychan or Buoy 7 = B7), using only physics-based abiotic covariates (light, salinity and hydrodynamics, as opposed to nutrient concentrations). 

Five different interaction matrices are considered: 

* ``null``/diagonal matrix
* ``unconstrained``/full matrix
* ``diatdin``, removing interactions between diatoms and dinoflagellates (interactions occur only within diatoms and within dinoflagellates)
* ``pencen``, also removing interactions between pennate and centric diatoms (interactions occur only within pennate diatoms, within centric diatoms and within dinoflagellates)
* ``inverse``, allowing interactions only between diatoms and dinoflagellates but neither within diatoms nor dinoflagellates 
