# generate RATTLE constraints for acetonitrile
# 
#

BEGIN{
    nstart=1
    nmol_start=1
    nmol=382
    natoms=6
    for(i = nmol_start; i <= nmol; i++){
        j=(i-1)*natoms + nstart
        printf "%8d %7d %14.7f %7d\n", j,   j+1, 1.090, 0
        printf "%8d %7d %14.7f %7d\n", j,   j+2, 1.090, 0
        printf "%8d %7d %14.7f %7d\n", j,   j+3, 1.090, 0
#       printf "%8d %7d %14.7f %7d\n", j,   j+4, 1.458, 0
        printf "%8d %7d %14.7f %7d\n", j+1, j+2, 1.78028, 0
        printf "%8d %7d %14.7f %7d\n", j+1, j+3, 1.78028, 0
        printf "%8d %7d %14.7f %7d\n", j+2, j+3, 1.78028, 0
#       printf "%8d %7d %14.7f %7d\n", j  , j+5, 2.61500, 0
#       printf "%8d %7d %14.7f %7d\n", j+4, j+5, 1.15700, 0
    }}
