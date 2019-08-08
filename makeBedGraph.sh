#!/bin/bash
folder=$1
path="/raidixshare_log-g/Dylan/190617_PE75/RepairSeq_062819/work/align/"
/home/preed/bin/makeTagDirectory $folder $path$folder"/"$folder"-sort.bam" -single
/home/preed/bin/makeUCSCfile $folder -o auto
