#!/bin/bash
for i in {2,4,8,16,32,64,128,256,512,1024,2048} 
do
	Rscript ~/qpMCMC/code/IsingExp.R -nProps $i &
done
