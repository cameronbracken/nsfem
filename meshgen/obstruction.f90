subroutine obstruction(h,k,radius,nholepts,nobspts,oshape,nholes,obstructxy)
    !!-------------------------------------------------------------------------
    !!obstruction Defines nholes holes in the problem domain
    !!
    !!
    !!Author:
    !!Cameron Bracken
    !!
    !!Date:
    !!19 Mar 2008
    !!
    !!~Outputs~
    !!
    !!obstructxy        -[double (2,nobspts)] the x,y coordinates of the
    !!                   circular obstructions
    !!
    !!~Inputs~
    !!
    !!h,k               - [double (nholes)] coordinates of the center of the 
    !!                    nholes circular holes in the domain
    !!radius            - [double (nholes)] the radius of each of the holes
    !!nholepts          - [integer] Number of holes/obstructions/circles
    !!nobspts           - [integer (nholes)] the number of points to use for 
    !!                    each of the holes
    !!oshape            - [character] the shape of the hole, anything other 
    !!                    than 'circle' is not currently implemented
    !!nholes            - [integer] the number of holes in the domain
    !!-------------------------------------------------------------------------`

    implicit none
    integer::nholes
    double precision,dimension(nholes)::h,k,radius
    double precision::thstep
    integer::i,j,hole,pcount,nobspts
    integer,dimension(nholes)::nholepts
    character(len=6)::oshape
    double precision,dimension(2,nobspts)::obstructxy
   

    if(oshape=="circle")then
        
        pcount=0
        do hole=1,nholes
            if(hole==1)then
                pcount=1
            else
                pcount=pcount+nholepts(hole-1)
            end if

            thstep=(2*3.1415926)/dble(nholepts(hole))
            
            j=0
            do i=pcount,pcount+nholepts(hole)-1
                j=j+1
                obstructxy(1,i)=radius(hole)*cos(dble(j)*thstep)+h(hole)
                obstructxy(2,i)=radius(hole)*sin(dble(j)*thstep)+k(hole)
            end do

        end do
        
    !elseif(oshape=="square")then
    !    
    !    obstructxy(1,1)=h-radius
    !    obstructxy(2,1)=k-radius
    !    obstructxy(1,2)=h-radius
    !    obstructxy(2,2)=k+radius
    !    obstructxy(1,3)=h+radius
    !    obstructxy(2,3)=k+radius
    !    obstructxy(1,4)=h+radius
    !    obstructxy(2,4)=k-radius
    !        
    end if 
    
    return
end
