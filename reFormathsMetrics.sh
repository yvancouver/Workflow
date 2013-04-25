cat $1| awk -F "\t" ' {
if ( NF != 0 ){
	FS Filed Separator
		if ( $1 ~/^#/) {
			print $0}
		else
			for (f = 1; f <= NF; f++)
				a[NR, f] = $f
}
else 
	print $0
}
NF > nf { nf = NF }
END {
	for (f = 1; f <= nf; f++)
	for (r = 1; r <= NR; r++)
		printf a[r, f] (r==NR ? RS : FS) 
}' | awk '{sub(/\t{6}/,"\t");print}'