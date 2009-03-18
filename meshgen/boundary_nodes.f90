subroutine boundary_nodes(boundary_xynodes,nsidepts,up,down,left,right)
    !!------------------------------------------------------------------------
    !!Generates nsidepts nodes along each side of a rectangle defined by 
    !!up,down,left,right. NOTE: A future step would be to generalize this to 
    !!allow for an arbitrarily shaped domain.
    !!
    !!Author:
    !!Cameron Bracken
    !!
    !!Date:
    !!19 Mar 2008
    !!
    !!~Outputs~
    !!
    !!boundary_xynodes      - [double (2,ninteriornodes)] x,y coordinates of 
    !!                        the nodes along the boundary
    !!
    !!~Inputs~
    !!
    !!nsidepts              - [integer] number of points along a side of the
    !!                        domain (the number of boundary points is) 
    !!                        4*nsidepts
    !!up,down,left,right    - [double] the edge coordinates of the boundary
    !!
    !!~Local Variables~
    !!
    !!xstep,ystep           -[double] distance between points along a 
    !!                       horizontal/vertical boundary 
    !!
    !!------------------------------------------------------------------------
    implicit none
    integer::nsidepts,i
    double precision,dimension(2,nsidepts*4)::boundary_xynodes
    double precision::up,down,left,right
    double precision::xstep,ystep
    
    boundary_xynodes(:,:)=0d0
    !define boundary points

    do i=1,nsidepts-1
        ystep=(up-down)/dble(nsidepts)
        xstep=(right-left)/dble(nsidepts)
        !up
        boundary_xynodes(1,i)=dble(i)*xstep
        boundary_xynodes(2,i)=up
        !down
        boundary_xynodes(1,nsidepts+i-1)=dble(i)*xstep
        boundary_xynodes(2,nsidepts+i-1)=down
        !left
        boundary_xynodes(1,2*nsidepts+i-2)=left
        boundary_xynodes(2,2*nsidepts+i-2)=dble(i)*ystep
        !right
        boundary_xynodes(1,3*nsidepts+i-3)=right
        boundary_xynodes(2,3*nsidepts+i-3)=dble(i)*ystep
    end do
        
        !Add the corners
        !upleft
        boundary_xynodes(1,4*nsidepts-3)=left
        boundary_xynodes(2,4*nsidepts-3)=up
        !downleft
        boundary_xynodes(1,4*nsidepts-2)=left
        boundary_xynodes(2,4*nsidepts-2)=down
        !downright
        boundary_xynodes(1,4*nsidepts-1)=right
        boundary_xynodes(2,4*nsidepts-1)=down
        !upright
        boundary_xynodes(1,4*nsidepts)=right
        boundary_xynodes(2,4*nsidepts)=up
    
    return
end
