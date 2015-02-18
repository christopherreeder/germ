#!/bin/bash

java -Xmx2g -cp /cluster/ccr/deconvolve/deconvolve.jar edu.mit.csail.cgs.reeder.parzen.BlindDeconvolve --numouterinters $numouterinters --numinnerinters $numinnerinters --finfile "$finfile" --cinfile "$cinfile" --ginfile "$ginfile" --clenfile "$clenfile" --numc $numc --goutfile "$goutfile" --foutfile "$foutfile" --coutfile "$coutfile"
