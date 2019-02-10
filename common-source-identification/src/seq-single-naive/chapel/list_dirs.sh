#!/bin/sh

for i in $(find $1 -maxdepth 1 -type d) ; do 
    echo -n $i": " ; 
    ( find $i -type f | wc -l ) ; 
done
