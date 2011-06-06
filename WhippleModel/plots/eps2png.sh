#! /bin/bash
for f in `ls *ps`; do
    convert -density 200 $f ${f%.*}.png;
done
