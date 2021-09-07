# connectivity matrix with generic central atom

BEGIN{
    if(natoms==""){
        print "ERROR: missing natoms"
        exit
    }
    if(nw==""){
        print "ERROR: missing nw"
        exit
    }
    if(nstart=="") nstart=1 
        for(i=nstart; i<=natoms; i=i+nw){
            printf "%6d",i
            for(j=1; j<nw; j++)
                printf "%6d %s",i+j, "1.0"
            printf "\n"
            for(j=1; j<nw; j++)
                printf "%6d\n",i+j
        }
}
