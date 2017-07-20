#!/bin/bash

make -j 4

mkdir bin

mv btrim bin

tar -czvf bin.tar.gz bin;


echo The end !;

