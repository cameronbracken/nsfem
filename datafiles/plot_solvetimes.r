slap=as.matrix(read.table('~/Desktop/codes/nsfem/datafiles/slaptime.out'))
lapack=as.matrix(read.table('~/Desktop/codes/nsfem/datafiles/lapacktime.out'))

pdf('~/Desktop/codes/nsfem/plots/solution_times.pdf',family='serif',pointsize=13)
plot(sort(lapack[,1]),lapack[order(lapack[,1]),2],type='l',col='red',xlab='# of Unknowns',ylab='Solution Time (s)',family='serif')
lines(sort(slap[,1]),slap[order(slap[,1]),2],col='blue')
legend('topleft',legend=c("Direct Dense Solver (LAPACK's DGETRS)","Iterative Sparse Solver (SLAP's DGMRES)"),lty='solid',col=c('red','blue'))
dev.off()

plot(sort(lapack[,1]),lapack[order(lapack[,1]),2],type='l',col='red',xlab='# of Unknowns',ylab='Solution Time (s)',family='serif')
lines(sort(slap[,1]),slap[order(slap[,1]),2],col='blue')
legend('topleft',legend=c("Direct Dense Solver (LAPACK's DGETRS)","Iterative Sparse Solver (SLAP's DGMRES)"),lty='solid',col=c('red','blue'))