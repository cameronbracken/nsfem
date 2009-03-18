nr=nrow(c)
nc=ncol(c)

quartz(width=10,height=10)

plot(1,xlab='Column',ylab='Row',xlim=c(0,nr),ylim=c(0,nc),type='n')

for(i in 1:nr){
	for(j in 1:nc){
		if(c[i,j]!=0){
			points(j,nr-i+1,cex=.2)
			}
		
		}
	
	}