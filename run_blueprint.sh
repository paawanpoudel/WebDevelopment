#!/bin/bash
#export PATH=$PATH:/usr/share/play_framework/pbin
cd /home/pp376/NetBeansProjects/blueprint

#R CMD Rserve --RS-conf /usr/share/blueprint/R_commands/Rserve.conf --RS-source /usr/share/blueprint/R_commands/plots_Progenitor_final.R

R CMD Rserve --RS-conf /home/pp376/RServe/rserve.conf   --RS-source /home/pp376/RServe/20140606_plots_Progenitor_final.R

nohup play run >| blueprint.log 2>| blueprint.err &
