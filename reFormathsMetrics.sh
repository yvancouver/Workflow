cat $1| awk -F "\t" ' { 													### NF number of field
if ( NF != 0 ) 										### RS 
{													### FS Filed Separator
if ( $1 ~/^#/) 										### test if line start with #
{
print $0
}
else												### 
for (f = 1; f <= NF; f++)					
a[NR, f] = $f										### NR
}
else 
print $0
}
NF > nf { nf = NF }
END {
for (f = 1; f <= nf; f++)
for (r = 1; r <= NR; r++)
printf a[r, f] (r==NR ? RS : FS) 
}'