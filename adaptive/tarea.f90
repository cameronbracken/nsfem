function tarea(xx,yy)
        implicit none
        double precision,dimension(3)::xx,yy
        double precision::tarea
        tarea=0.5d0*dabs(xx(1)*(yy(2)-yy(3)) + xx(2)*(yy(3)-yy(1)) + xx(3)*(yy(1)-yy(2)))
end