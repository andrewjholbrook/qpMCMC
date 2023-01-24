#!/bin/bash
for i in {1..100} 
do
	Rscript ~/qpMCMC/code/pMCMCvsQPMCMC.R -seed $i &
done
