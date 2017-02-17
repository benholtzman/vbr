#!/bin/bash

mkdir source

cp -r 01* source/
cp -r 02* source/
cp -r 03* source/

mkdir source/04_scripts
cp -r 04_scripts/*.sh source/04_scripts/
cp TwoPhase*.m source/
cp $2 source/

tar -cvzf source.tar.gz source/ > /dev/null

rm -r source
mv source.tar.gz ${1}/
