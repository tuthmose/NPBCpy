# Generate fragment string for G16 input
#

BEGIN{
	if(nfrag==""){
		print "ERROR: missing nfrag"
		stop=1
		exit
	}
	if(natoms==""){
		print "ERROR: missing natoms"
		stop=1
		exit
	}
       
        if(nstart=="") nstart=1
        for(i=nstart;i<=nfrag;i++){
        	for(j=1;j<=natoms;j++){
			printf "\(fragment=%5d\)\n",-i
		}
	}}
