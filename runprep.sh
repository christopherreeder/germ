#!/bin/bash

java -Xmx4g -cp /cluster/ccr/deconvolve/deconvolve.jar edu.mit.csail.cgs.reeder.parzen.BlindDeconvolvePrep --genomefile "$genomefile" --targetsize $targetsize --h $h --binsize $binsize --readfile "$readfile" --region "$region" --selfliglikfile "$selfliglikfile" --storagesize $storagesize --foutbase "$foutbase" --coutbase "$coutbase" --clenoutbase "$clenoutbase"
