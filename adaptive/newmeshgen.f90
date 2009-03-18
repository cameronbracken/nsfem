subroutine newmeshgen(xynodes,ntotnodes,file,filename,verbose)
    !!-----------------------------------------------------------------------       
    !!The purpose of this subroutine is to write a file containing node 
    !!information to a .poly file that is accepted by 'triangle'.
    !! 
    !!This is almost identical to the initial mesh generation 
    !!
    !!Local Variables:
    !!  
    !!nnodes        - [integer] total number of nodes in domain
    !!xynodes       - [double (nnodes,2)]  x,y coordinates all of nodes
    !!                includes boundary nodes interior nodes and hole nodes
    !!obstructxy    - [double (2,ninteriornodes)] x,y coordinates of the 
    !!               nodes around the hole 
    !!theend        - [integer] the number of node coordinates to print out 
    !!              if verbose=.true.
    !!-----------------------------------------------------------------------
    implicit none
    double precision::up,down,left,right
    double precision,dimension(ntotnodes,2)::xynodes
    integer::nnodes,i,j,nobspts,ntotnodes,nholes
    integer::nsidepts,ninteriornodes,hole,pcount
    integer,allocatable,dimension(:)::index,nholepts
    double precision,dimension(:),allocatable::holeh,holek,holer
    character(len=*)::filename,file
    character(len=6)::oshape
    logical::verbose
    
    open(15,file=file)
    read(15,*)up,down,left,right,ninteriornodes,nsidepts,nobspts,nholes
    close(15)   

    open(15,file='holefile')
    allocate(holeh(nholes),holek(nholes),holer(nholes),nholepts(nholes))
    do i=1,nholes
        read(15,*)holeh(i),holek(i),holer(i),nholepts(i),oshape
    end do
    close(15)
                                            

    if(verbose)then
        write(*,*)' '
        write(*,*)'There are',ntotnodes,'nodes'
        write(*,*)' '
        write(*,*)'node coordinates'
        do i=1,ntotnodes
            write(*,'(1i5,2f10.3)')i,xynodes(i,:)
        end do
        write(*,*)'     .         .   '
        write(*,*)'     .         .   '
        write(*,*)'     .         .   '
        write(*,*)ntotnodes-ntotnodes,'more'
        write(*,*)' '
    endif
    
    
    
    !--------make .poly file------------------------------
    open(12,file=filename)
    write(12,*)ntotnodes,2,0
        !node section header
        !number of nodes,dimesion,specify boundary
    do i=1,ntotnodes
            write(12,"(1i6,2f25.15,3i4)")i,xynodes(i,1),xynodes(i,2)
    end do
    
    if(oshape/='none')then
        write(12,*)nobspts,0
            !polygon section header
            !number of polygon points,specify boundary
        
        do hole=1,nholes
	        if(hole==1)then
	            pcount=1
	        else
	            pcount=pcount+nholepts(hole-1)
	        end if

	        do i=pcount,pcount+nholepts(hole)-1
	            if(i<(pcount+nholepts(hole)-1))then
	                    !all but one of the polygon lines
                    write(12,*)i,ninteriornodes+4*nsidepts+i,ninteriornodes+4*nsidepts+i+1
	            else
	                    !need to connect the last point to the first one on the last step
                    write(12,*)i,ninteriornodes+4*nsidepts+i,ninteriornodes+4*nsidepts+pcount
	            end if
	        end do
	    end do
        
        write(12,*)nholes
            !hole header
        do i=1,nholes
            write(12,*)i,holeh(i),holek(i)
        end do
            !hole number, x,y interior point
    endif
    
    close(11)
    close(12)

    if(verbose)then
        write(*,*)' '
        write(*,*)"wrote node and polygon (obstruction) information to ",filename
        write(*,*)' '
    endif
        
end
