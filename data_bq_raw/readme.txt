Raw images have been submitted to XXXXXXX

#### Data structure ####
Each directory contains the data from a plate in three seperate replicates (runs).
The data contains the calculated global biofilm parameters from BiofilmQ as well as structural information about
the wells and experiment.

-	Well corresponds to the microwell of the titerplate (A1 is the first well, top left).
-	Well_pos corresponds to a focal position in the XY-direction in a single well 
	(1-5, as there are five measurement positions). This is equivalent to a technical replicate.
-	Ch corresponds to the channel. ch = 1 is red, ch = 2 is green.
-	Run corresponds to the biological replicate number. Each run was taken on a seperate day.

IN ITS CURRENT STATE, RUN AND WELL_POS ARE SWITCHED!