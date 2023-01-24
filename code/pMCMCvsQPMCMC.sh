#!/bin/bash
for i in {1..50} 
do
	Rscript ~/qpMCMC/code/pMCMCvsQPMCMC.R -seed $i &
done
