# G Mancini April 2021
# 
# calculate GLOB U_{vdW} potential
# from two column file with layer number and 
# gaussian number of each optmization cycle

BEGIN{
	if(last=="")
		last="last.dat"
	if(cumulative=="")
		cumulative="cumulative.dat"
        nlayer = 81
	nskip= 10
	sigma  = 0.125
	sigma2 = 2.0*sigma*sigma
	height = 0.010
	npoints = 200
	Rmax = 20.25
        step = 0.1
}

nskip > 0 && nread <= nskip{
	layer_i = NR % nlayer
        if(layer_i == 1)
        	nread++
}

nskip <= 0 || (nread > nskip){
	layer_i = NR % nlayer
        if(layer_i == 1) nmeta++
        if(layer_i == 0) layer_i = nlayer
        tmp[layer_i] = $NF - tmp[1]
        tmp_acc[layer_i] += tmp[layer_i]
        radh[layer_i] = $(NF-1)
}

END{
	print nmeta," optmization cycles read"
	for(i=1;i<=nlayer;i++)
            tmp_acc[i] = tmp_acc[i]/nmeta
	for(i=npoints;i>=1;i--)
        {
            Sum = 0.0
            Sum_last = 0.0 
            RO = i*step
	    for(j=1;j<=nlayer;j++)
            {
		R2 = (RO - radh[j])**2
		earg = R2/sigma2
                if(earg < 100.0)
                {
                        gd = (tmp_acc[j] - tmp_acc[1])
                        H  = height*gd
			Sum += H*exp(-earg)
                        gd = (tmp[j] - tmp[1])
                        H  = height*gd
			Sum_last += H*exp(-earg)

                }}
       	    printf "%12.6f %12.6f\n",Rmax-RO,Sum > cumulative
       	    printf "%12.6f %12.6f\n",Rmax-RO,Sum_last > last
}}
