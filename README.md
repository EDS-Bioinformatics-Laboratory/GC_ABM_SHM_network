# Publication
_In preparation_

This software accompanies the publication

**Understanding repertoire sequencing data through a computational model of the germinal center.**

Rodrigo García-Valiente, Elena Merino Tejero, Maria Stratigopoulou, Daria Balashova, Aldo Jongejan, Danial Lashgari, Aurélien Pélissier, Tom G. Caniels,
Mathieu A. F. Claireaux, Aram Al-Soudi, Marit J. van Gils, María Rodríguez Martínez, Niek de Vries, Michael Meyer-Hermann, Jeroen E.J. Guikema, Huub Hoefsloot,Antoine H.C. van Kampen

# Project
## GC_ABM_SHM_network

This repository includes code ONLY of the multiscale model of plasma cell diferentiation + sequence representation + SHM fate tree in germinal centers. The Agent-based model is based on Mafalda (Merino-Tejero, E., 2021), which is based on [Hyphasma](https://www.helmholtz-hzi.de/en/research/research-topics/immune-response/systems-immunology/our-research/) (e.g., Michael-Meyer Hermann, 2012).  
The GRN is based on Martinez et al., 2012. 

## How-to run
-Modify base_Path variable.
-Modify /home/rgarcia/Escritorio/NGly_scripts/Control_VJ_GC_model path to the path where your NT Fab multifasta (separated per region) files are.
-Create a Sequence_db.csv in bcinflow09 where your whole Fab input AA sequences are assigned a coordinate in the shape space (see example in folder).
-Give it a try!

## Software
All software is written in C++

## References
* Merino-Tejero, E., Lashgari, D., García-Valiente, R., Gao, X., Crauste, F., Robert, P. A., Meyer-Hermann, M., Rodríguez-Martínez, M. , van Ham, S. M., Guikema, J. E. J., Hoefsloot, H., van Kampen, A. H. C. (2021). Multiscale Modeling of Germinal Center Recapitulates the Temporal Transition From Memory B Cells to Plasma Cells Differentiation as Regulated by Antigen Affinity-Based Tfh Cell Help. Frontiers in Immunology, 11.
* Martínez, M. R., Corradin, A., Klein, U., Álvarez, M. J., Toffolo, G. M., di Camillo, B., … Stolovitzky, G. A. (2012). Quantitative modeling of the terminal differentiation of B-cells and mechanisms of lymphomagenesis. Proceedings of the National Academy of Sciences, 109(7), 2672–2677. 
* Meyer-Hermann, M., Mohr, E., Pelletier, N., Zhang, Y., Victora, G. D., & Toellner, K. M. (2012). A theory of germinal center b cell selection, division, and exit. Cell Reports, 2(1), 162–174. 
