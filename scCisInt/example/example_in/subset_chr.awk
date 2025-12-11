awk -v chr="$1" '
NR==1 {
    keep[1]=1
    cols[1]=1
    n=1
    for(i=2;i<=NF;i++){
        if($i ~ chr){
            n++
            cols[n]=i
        }
    }
    for(j=1;j<=n;j++)
        printf "%s%s", $cols[j], (j<n?OFS:ORS)
    next
}
{
    for(j=1;j<=n;j++)
        printf "%s%s", $cols[j], (j<n?OFS:ORS)
}
' "$2"

