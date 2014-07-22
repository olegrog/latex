#!/bin/bash

conv() {
	printf "Converting $f..."
	file=$(echo $f | sed 's/\.pdf//')
	pdfcrop --margins=20 $f $f > /dev/null 2>&1
	[ -f $file.png ] || gs -sOutputFile=$file.png -sDEVICE='pngalpha' -r300 -dBATCH -dNOPAUSE $file.pdf
	convert -transparent white -resize 120x120 $file.png $file.png
	echo "done"
}

[[ $# -ne 0 ]] && list="$@" || list=$(ls | grep pdf)

for f in $list; do
	conv
done
