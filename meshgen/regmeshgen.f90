program meshgen
    !!-----------------------------------------------------------------------    
    !!The purpose of this program is to write a file containing node 
    !!information to a .poly  file that is accepted by the triangle meshing
    !!software.  This subroutine generates one set of nodes based on the 
    !!parameters of the domain.  A file 'domain' is written describing the domain
    !!A few subroutines are called.  See those for additioanl documentation 
    !!
    !!~Domain Parameters(in order of appearance)~
    !!NOTE: These parameters are set up for a rectangular domain
    !!
    !!oshape            - [character] the shape of the obstruction in the domain
    !!                    values include: 'circle','square','none'
    !!up                - [double] y value of the top boundary
    !!down              - [double] y value of the bottom boundary
    !!left              - [double] x value of the left boundary
    !!right             - [double] x value of the right boundary
    !!h,k               - [double] centroid point of the obstruction
    !!radius            - [double]  radius of the obstruction , for 'square' 
    !!                    radius is the radius of the inscribed circle
    !!filename          - [character] Name of the output file 
    !!                                  (command line input)
    !!
    !!~Mesh Parameters(in order of appearance)~
    !!
    !!nsidepts          - [integer] number of points along an outside boundary.
    !!                    NEEDS TO BE GENERALIZED FOR RECTANGLES
    !!nobspts           - [integer] The number of pts used to define the 
    !!                    obstruction. Only makes sense for 'circle'
    !!                        oshape='square' ----> nobspts=4
    !!                        oshape='none'   ----> nobspts=0
    !!ninteriornodes    - [integer] NEEDS TO BE CHANGED
    !!
    !!~Local Variables~
    !!
    !!verbose           - [logical] should hella crap be printed out
    !!nnodes            - [integer] total number of nodes in domain
    !!nrows             - [integer] number of regular mesh columns 
    !!ncols             - [integer] number of regular mesh rows
    !!xynodes           - [double (2,nnodes)]  x,y coordinates all of nodes
    !!                    includes boundary nodes interior nodes and hole nodes
    !!                    note: x-----> 
    !!                          y----->  by rows
    !!interior_xynodes  - [double (2,ninteriornodes)] x,y coordinates of 
    !!                    the interior nodes, if any 
    !!boundary_xynodes  - [double (2,ninteriornodes)] x,y coordinates of 
    !!                    the nodes along the boundary
    !!obstructxy        - [double (2,ninteriornodes)] x,y coordinates of the 
    !!                    nodes around the hole 
    !!-----------------------------------------------------------------------
    implicit none
    double precision::up,down,left,right,radius,h,k,rand,lh
    double precision,allocatable,dimension(:,:)::interior_xynodes,&
                                        boundary_xynodes,xynodes,obstructxy
    integer::nnodes,i,j,theend,pcount,hole
    integer::nsidepts,ninteriornodes,nrows,ncols,nignore,nholes,nobspts
    integer,allocatable,dimension(:)::index,nholepts
    double precision,dimension(:),allocatable::holeh,holek,holer
    character(len=20)::filename
    character(len=6)::oshape
    logical::verbose,holes
    
    call getarg(1,filename)
   
    open(25,file='domain')
    read(25,*)up,down,left,right,nsidepts,nrows,ncols,verbose,holes
    close(25)

    nobspts=0
    if(holes)then
        open(26,file='holefile')
        read(26,*)nholes
        write(*,*)'there are ',nholes,'holes'
        allocate(holeh(nholes),holek(nholes),holer(nholes),nholepts(nholes))
        do i=1,nholes
            read(26,*)holeh(i),holek(i),holer(i),nholepts(i),oshape
            nobspts=nobspts+nholepts(i)
        end do
    end if 
    
    
    allocate(obstructxy(2,nobspts))
    
    call obstruction(holeh,holek,holer,nholepts,nobspts,oshape,nholes,obstructxy)
    
    
    allocate(interior_xynodes(2,nrows*ncols))


        !generates a set of interior points
    call interior_nodes(nignore,interior_xynodes,nrows,ncols,up,down,&
                            left,right,nholes,h,k,radius)
    ninteriornodes=nrows*ncols-nignore
    
    nnodes=(ninteriornodes+4*nsidepts+nobspts)

    allocate(boundary_xynodes(2,4*nsidepts),xynodes(2,nnodes))
    
        !generates the boundary and corner nodes
    call boundary_nodes(boundary_xynodes,nsidepts,up,down,left,right)
        
    !write(*,*)
    !do i=1,4*nsidepts
    !   write(*,'(2f10.3)')boundary_xynodes(1,i),boundary_xynodes(2,i)
    !end do
    
        !stick the boundary and interior nodes together 
    xynodes(:,1:ninteriornodes)=interior_xynodes
    xynodes(:,(ninteriornodes+1):(ninteriornodes+4*nsidepts))=boundary_xynodes
    xynodes(:,(nnodes-nobspts+1):nnodes)=obstructxy
    
    theend=nnodes

    if(verbose)then
        write(*,*)' '
        write(*,*)'There are',nnodes,'nodes'
        write(*,*)' '
        write(*,*)'node coordinates'
        do i=1,theend
            write(*,'(2f10.3)')xynodes(1,i),xynodes(2,i)
        end do
        write(*,*)'     .         .   '
        write(*,*)'     .         .   '
        write(*,*)'     .         .   '
        write(*,*)nnodes-theend,'more'
        write(*,*)' '
    endif
    
    
    !--------make .poly file------------------------------
    open(12,file=filename)
    write(12,*)nnodes,2,0
        !node section header
        !number of nodes,dimesion,specify boundary
    do i=1,nnodes
        write(12,"(1i3,2f10.4,3i4)")i,(xynodes(j,i),j=1,2)
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
                    write(12,*)i,nnodes-nobspts+i,nnodes-nobspts+i+1
                        !all but one of the polygon lines
                else
                        !need to connect the last point to the first one on the last step
                    write(12,*)i,nnodes-nobspts+i,nnodes-nobspts+pcount
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
        !write(*,*)"wrote node information to ",nodefile    
        !write(*,*)' '
        write(*,*)' '
    endif
        
    !    write(*,*)"wrote node and polygon (obstruction) information to ",filename

        open(15,file='domain.out')
        write(15,*)up,down,left,right,ninteriornodes,nsidepts,nobspts,nholes
        close(15)
end
