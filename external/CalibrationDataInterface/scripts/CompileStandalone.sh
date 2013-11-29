#!/bin/sh

for dir in obj obj/dep ; do
  if ! [ -d $dir ] ; then
    mkdir $dir
  fi
done

#ln -s ../CalibrationDataInterface .
for file in CalibrationDataInterfaceBase.cxx CalibrationDataInterfaceROOT.cxx CalibrationDataContainer.cxx CalibrationDataEigenVariations.cxx ; do
  sed -e "s/CalibrationDataInterface\///g" ../Root/$file > ./$file
done
for file in CalibrationDataInterfaceBase.h CalibrationDataInterfaceROOT.h CalibrationDataContainer.h CalibrationDataVariables.h CalibrationDataEigenVariations.h ; do
  sed -e "s/CalibrationDataInterface\///g" ../CalibrationDataInterface/$file > ./$file
done

make
