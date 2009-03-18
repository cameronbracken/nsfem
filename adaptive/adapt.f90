program element_error
    implicit none 
    character(len=100)::node_file,ele_file,solution_file,corner_file,outputfile,error_method
    double precision,dimension(:,:),allocatable::node_xy,indicator,newnodes,allnodes
    integer,dimension(:,:),allocatable::ele_nums
    integer::err,linecount,i,j,k,nnode,nele,ncol
    integer::point,nrefine
    double precision::x,y,p,minv
    double precision,dimension(:),allocatable::ele_err,z,line
    double precision,dimension(3)::xx,yy
    integer,dimension(:),allocatable::refine,corner
    double precision::tarea,area,globalmax,error
    integer::ncorner,ern
    logical::verbose,quiet
    double precision::minarea,maxrelerr,toler
    integer::maxviolations
    !logical,parameter::verbose=.false.,quiet=.true.
    !double precision,parameter::minarea=0.001d0,maxrelerr=0.15d0,toler=0.001d0
    !integer,parameter::maxviolations=3

    open(26,file='refinefile')
    read(26,*)error_method,verbose,quiet,minarea,maxrelerr,toler,maxviolations,ern
    !write(*,*)verbose,quiet,minarea,maxrelerr,toler,maxviolations,ern
    close(26)

    !!this program uses the output from a finite element model
    !!along with the element and node information to estimate 
    !!elemental error
   
    write(*,*)
    write(*,*)'Element error analysis...'   
    
    call getarg(1,node_file)
    call getarg(2,ele_file)
    call getarg(3,solution_file)
    call getarg(4,corner_file)
    call getarg(5,outputfile)
    
    !--------------------------------------------------------------------
    if(.not.quiet)then 
        write(*,*)'Reading ',trim(corner_file),'...'
    end if
    linecount=0
    call itable_header_read(corner_file,ncol,linecount)
    allocate(corner(linecount))

    open(14,file=corner_file)
    ncorner=0
    do i=1,linecount
        read(14,*)corner(i)
        if(corner(i)==1)then
            ncorner=ncorner+1
        endif
    end do
    close(14)
    nnode=ncorner
    if(.not.quiet)then 
        write(*,*)'Read',linecount,'lines and there were',ncorner,'vertex nodes'
        write(*,*)
    end if
    !call print_double_mat(corner,linecount,1)


    !--------------------------------------------------------------------
    if(.not.quiet)then 
        write(*,*)'Reading ',trim(node_file),'...'
    end if
    open(11,file=node_file)

    allocate(node_xy(ncorner,2))
    j=0
    do i=1,linecount
        read(11,*)x,y
        if(corner(i)==1)then
            j=j+1
            node_xy(j,1)=x
            node_xy(j,2)=y
        endif
    end do
    close(11)
    if(.not.quiet)then 
        write(*,*)'Read',ncorner,'lines out of',linecount
        write(*,*)
    end if
    !call print_double_mat(node_xy,ncorner,2)

    
    !--------------------------------------------------------------------
    if(.not.quiet)then 
        write(*,*)'Reading ',trim(solution_file),'...'
    end if
    
    
    call dtable_header_read(solution_file,ncol,linecount)
    !write(*,*)ncol,linecount
    !allocate(indicator(ncorner,ncol))                                                                           
    !call read_dtable(solution_file,linecount,ncol,0,indicator) 

    open(13,file=solution_file)
    allocate(line(ncol),indicator(ncorner,3))
    j=0
    do i=1,linecount
        read(13,*)line(:)
       ! write(*,*)line
        if(corner(i)==1)then
            j=j+1
            indicator(j,:)=line(:)
        endif
    end do
    close(13)
    
    if(.not.quiet)then 
        write(*,*)'Read',ncorner,'lines out of',linecount
        write(*,*)
    end if
    !call print_double_mat(indicator,ncorner,2)
    
    
    !--------------------------------------------------------------------
    if(.not.quiet)then 
        write(*,*)'Reading ',trim(ele_file),'...'
    end if
    linecount=0
    call itable_header_read(ele_file,ncol,linecount)
    allocate(ele_nums(linecount,ncol))
    call read_itable(ele_file,linecount,ncol,0,ele_nums)
    
    nele=linecount
    if(.not.quiet)then 
        write(*,*)'Read',linecount,'lines'
        write(*,*)
    end if
    !call print_int_mat(ele_nums,linecount,3)

    !--------------------------------------------------------------------

    allocate(ele_err(nele),refine(nele),z(3))
    refine=0
    nrefine=0
    
    if(.not.quiet)then
        write(*,*)'The error indicator is ',trim(error_method)
        write(*,*)  
    end if

    globalmax=0d0
    do i=1,nele
         
        z(:)=indicator(ele_nums(i,:),ern)
        xx(:)=node_xy(ele_nums(i,:),1)
        yy(:)=node_xy(ele_nums(i,:),2)  
        
        if(dabs(maxval(z)-minval(z))>globalmax)globalmax=dabs(maxval(z)-minval(z))

    end do        

    !globalmax=maxval(dabs(indicator(:,ern)))
    !globalmax=300
    write(*,*)'maximum indicator value is',globalmax

    do i=1,nele
    
        !get the node coords and solution values for this element
        do k=1,3
            point=ele_nums(i,k)
            z(k)=indicator(point,ern)
            xx(k)=node_xy(ele_nums(i,k),1)
            yy(k)=node_xy(ele_nums(i,k),2)
        end do

        ele_err(i)=error(xx,yy,z,globalmax,toler,error_method)

        !write(*,*)'error',i,ele_err(i)
        
        area=tarea(xx,yy)   
        
        !write(*,*)'area',i,area

        if(ele_err(i)>maxrelerr .and. area>minarea)then
            refine(i)=1
            nrefine=nrefine+1
        end if
    
    end do
    
    write(*,*)'Average relative error is',sum(ele_err)/dble(nele) 
    write(*,*)'Minimum relative error is',minval(ele_err)  
    write(*,*)'Need to refine',nrefine,'elements'
    open(22,file='summary.out')
    write(22,*)'Average relative error is',sum(ele_err)/dble(nele) 
    write(22,*)'Minimum relative error is',minval(ele_err)  
    write(22,*)'Need to refine',nrefine,'elements'
    close(22)
    

        
    if(nrefine>maxviolations)then
        allocate(newnodes(nrefine,2),allnodes(nrefine+nnode,2))
        
            !refine the specified elements by placing a new node at their centre    
        call getnewnodes(refine,nrefine,'domain',node_xy,nnode,ele_nums,nele,allnodes,newnodes)
            !write the new mwsh file 
        call newmeshgen(allnodes,nnode+nrefine,'domain',outputfile,.false.)

        if(verbose)then
            open(14,file='error.out')
            do i=1,nele
                    write(14,*)ele_err(i)
            end do
            close(14)
            write(*,*)'Wrote elemental error to error.out'
        end if
        if(.not.quiet)then
            write(*,*)'Finished assembleing coeficient matrix.'
            write(*,*)
        end if
    
    else

        open(15,file='continue',status="old")
        write(15,*)'done'

    end if

end program element_error
