fn = 'rel_error_1320.out'

u = scan(fn,nlines=1)
v = scan(fn,nlines=1,skip=1)
p = scan(fn,nlines=1,skip=2)

n = max(length(u),length(v),length(p))

pdf('../plots/errorstep.pdf',family='serif',pointsize=13)
plot(u*100,col=2,xlab='Mesh Refinement Step',ylab='Average Percent Realative Elemental Error',ylim=c(100*range(u,v,p)),xlim=c(1,n),type='b')
lines(v*100,col=3,type='b')
lines(p*100,col=4,type='b')
name = c('Horizontal velocity','Vertical velocity','Pressure')
legend('topright',name,col=2:4,lty='solid')
dev.off()

plot(u*100,col=2,xlab='Mesh Refinement Step',ylab='Average Percent Realative Elemental Error',ylim=c(100*range(u,v,p)),xlim=c(1,n),type='b')
lines(v*100,col=3,type='b')
lines(p*100,col=4,type='b')
name = c('Horizontal velocity','Vertical velocity','Pressure')
legend('topright',name,col=2:4,lty='solid')