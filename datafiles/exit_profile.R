parms = scan('../domain',nlines=1,what='character')
up = as.numeric(substr(parms,1,1))

sol = as.matrix(read.table('solution.out'))
u = sol[,1]
v = sol[,2]
mag = sqrt(u^2+v^2)

xy = as.matrix(read.table('nodes.out'))

outmag = mag[xy[,2]==up]
print(outmag)
inmag = mag[xy[,2]==0]
xy = xy[xy[,2]==up,]
x = xy[,1]

pdf('../plots/vel_profile.pdf',family='serif',pointsize=13)
plot(sort(x),inmag[order(x)],type='b',col='steelblue',xlab='X',ylab='Velocity Magnitude')
lines(sort(x),outmag[order(x)],type='b')

legend('topright',c('Inlet','Outlet'),lty='solid',col=c('steelblue','black'))
dev.off()




plot(sort(x),inmag[order(x)],type='b',col='steelblue',xlab='X',ylab='Velocity Magnitude')
lines(sort(x),outmag[order(x)],type='b')

legend('topright',c('Inlet','Outlet'),lty='solid',col=c('steelblue','black'))