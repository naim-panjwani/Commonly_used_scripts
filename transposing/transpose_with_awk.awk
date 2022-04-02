gawk '
{ for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' /home/Bowei/tryImputation/chr1/concord/400dataALLSNPsameOrderwithRef.txt >bowei_transpose_with_awk.txt
