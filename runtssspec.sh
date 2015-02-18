#!/bin/bash

java -Xmx10g -cp /cluster/ccr/germ.jar edu.mit.csail.cgs.reeder.germ.InteractionSpecificity --genomefile "$genomefile" --pointfile "$pointfile" --tssfile "$tssfile" --binsize $binsize --margfile "$margfile" --interthresh $interthresh --backval $backval --numsamps $numsamps --readdistfile "$readdistfile" --readfile "$readfile" --storagesize $storagesize --outfile "$outfile" --runnum $runnum --perrun $perrun --jointsumfile "$jointsumfile"
