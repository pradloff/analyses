#!/bin/sh

for i in `ls *Check*.eps` ; do
  gv $i
done
