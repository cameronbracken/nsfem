function error(xx,yy,z,globalmax,toler,error_method)
        implicit none
        double precision,dimension(3)::xx,yy,z
        double precision::globalmax,toler,error
        character(len=*)::error_method
        
        if(trim(error_method)=='heuristic')then 
            if(maxval(dabs(z))<toler)then
                error=0d0
            else
                error=dabs(maxval(z)-minval(z))/globalmax
            endif
        end if
end
