#!/bin/bash

for f in $(ls | grep pdf)
do
	pdfcrop -margins 5 $f $f
done

