#!/bin/bash
for i in 1 2 3 4 5
do
    cat PhOptFullROCSMscan | ../bin/pXar -d ../data/testPhR0219/ -v DEBUGAPI
    cp ../data/testPhR0219/pxar.root ../data/testPhR0219/pxarWR$i.root 
done