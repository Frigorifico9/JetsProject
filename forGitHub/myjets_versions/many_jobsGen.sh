#!/bin/bash

for ((INDEX = 0; INDEX < 400; INDEX++)) 
do
   cd /storage/alice/frigorifico/jetsProject/allSettings
   rm seetings$INDEX.cmnd
   for i in seetings100000.cmnd  ######el archivo base
   do
   sed 's/thisGetsChanged/'$INDEX'/' seetings100000.cmnd >> seetings$INDEX.cmnd 
   done
cd /storage/alice/frigorifico/jetsProject/forCluster
   qsub -v JOBID=$INDEX generaEventosc.pbs
#   sleep 1  
done

