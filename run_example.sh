#!/bin/bash

chmod +x EasySpeciesTree.py
cd example/
for i in `ls *.fas`;do cat $i >>all.pep.fas;done
../EasySpeciesTree.py -in1 species_id.txt -in2 SingleCopyOrthogroups.txt -in3 Orthogroups.csv -in4 all.pep.fas -nb 10 -t 12
cd ../
