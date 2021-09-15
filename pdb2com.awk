# G Mancini mar 2021

## Assemble a PDB (with the right coordinates) and additional files in a G16 input

## natoms: number of atoms in charges and connectivity matrix
## cent:   box dimension / 2
## nstart: start printing connectivy and charges from this atom onward
## rigid:  print >0 or <0 propagation flag
## ARGIND=1 atoms and charges
## ARGIND=2 pdb
## ARGIND=3 connectivity matrix

BEGIN{
    if(natoms=="" || cent==""){
        print "ERROR: missing natoms in residue and center"
        stop=1
        exit
    }
    fmt1="%25-s%8d%16.6f%16.6f%16.6f"
    fmt2="%25-s%8d%16.6f%16.6f%16.6f\n"
    nmol=0
    atom_counter=0
    if(nstart == "")
        nmol = 0
    else if(nstart > 1)
        nmol = 1
}
        
ARGIND==1{
    A[FNR] = $1
    }
        
ARGIND==2 && $0~/ATOM/ {
    atom_counter++
    if (NF==11 && $NF=="0.00"){
        atx = 7
        aty = 8
        atz = 9
    } else if(NF==10 && $NF=="0.00"){
        atx = 6
        aty = 7
        atz = 8
    } else if(NF==11 && $NF~/[A-Z]/){
        atx = 6
        aty = 7
        atz = 8
    } else if(NF==10 && $NF~/[A-Z]/){
        atx = 5
        aty = 6
        atz = 7
    }
    if (atom_counter >= nstart){
        at = (atom_counter -nstart +1) % natoms
        if(at==1)
            nmol++
        else if(at==0)
            at=natoms
        if(rigid != "")
            snmol = -nmol
        else
            snmol = nmol
        frg="(fragment=" nmol ")"
        atom=A[at]
        a_q_frg = atom frg
        record[atom_counter] = sprintf(fmt1,a_q_frg,snmol,$atx - cent,$aty - cent,$atz - cent)
    } else {
        printf fmt2,$3,-1,$atx - cent,$aty - cent,$atz - cent
       }
}

ARGIND==3{
    row[FNR]=$0
}

END{
    if(stop)
        exit
    for(r=nstart; r<=atom_counter; r++)
        print record[r]
    printf "\n"
    for(i=1; i<nmol; i++)
        for(j=1; j<=natoms; j++){
            if(length(row[j] > 1)){
                a=split(row[j], mat )
                printf "%d  ",(i-1)*natoms+mat[1]+nstart-1
                for(k=2; k<a; k=k+2)
                    printf "%d  %2.1f  ",(i-1)*natoms+mat[k]+nstart-1,mat[k+1]
            printf "\n"
        }else{
            print row[j]
     }}
    printf "\n"
}
