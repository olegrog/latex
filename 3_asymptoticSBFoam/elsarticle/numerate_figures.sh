#!/bin/bash

[[ $# -ne 1 ]] && { echo "Usage: `basename $0` file.tex"; exit 1; }
[[ -f $1 ]] || { echo "No such file $1"; exit 1; }

old=".$1.bak"
new="$1"
cp $new $old

awk '
BEGIN { n = 1; }

$0 ~ /includegraphics/ {
	if (match($0, "{.*}")) {
		before = substr($0,1,RSTART);
		pattern = substr($0,RSTART+1,RLENGTH-2);
		after = substr($0,RSTART+RLENGTH-1);
		new = "Fig"n;
		printf("%s%s%s\n", before, new, after);
		system("mv "pattern".pdf "new".pdf");
		n++;
	}
}
$0 !~ /includegraphics/ { 
	print 
}
' $old  > $new
