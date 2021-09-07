# G Mancini Apr 2021
# get GLOB data and print it in a two columns file

BEGIN{
	if(nlayer==""){
		print "Error: set number of layer with awk -v nlayer=<x> ..."
		exit
	}
        if(last=="")
		last = 0
	thick = 0.25
}

/radius, group, density \(gr\.\), \#gau/{
	for(i=1; i<=nlayer; i++){
		getline
                if(NF==5)
                    print $2,$NF
                else if(NF==6)
                    print $3,$NF
                else{
                    print "ERROR parsing GLOB updates"
                    exit -1
	}}
}

