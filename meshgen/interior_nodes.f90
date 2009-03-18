subroutine interior_nodes(nignore,xynodes,nrows,ncols,up,down,left,right,nholes,h,k,r)
    !!------------------------------------------------------------------------
    !!Generates a uniform grid of points such that none of the points fall 
    !!inside of nholes circular shapes with in the grid.
    !!
    !!Author:
    !!Cameron Bracken 
    !!
    !!Date
    !!19 Mar 2008
    !!
    !!~Inputs~
    !!
    !!nrows,ncols       - [integer] number of rows/colums of gridpoints 
    !!up,down,left,right- [double] the edge coordinates of the boundary
    !!nholes            - [integer] number of hles in the domain
    !!h,k               - [double (nholes)] coordinates of the center of the 
    !!                    nholes circular holes in the domain
    !!r                 - [double (nholes)] the radius of each of the holes
    !!
    !!~Outputs~
    !!
    !!nignore           - [integer] the number of grid points which get ignored
    !!xynodes           - [double (2,nrows*ncols)] the x,y coordinates of
    !!                    the interior nodes which are created
    !!
    !!~Local Variables~
    !!
    !!rinc,cinc         - [double] the distance between rows/columns of 
    !!                    interior grid points
    implicit none
    integer::nholes
    double precision,dimension(2,nrows*ncols)::xynodes
    double precision::up,down,left,right,x,y,lh,rinc,cinc
    double precision,dimension(nholes)::h,k,r
    integer::nnodes,i,j,ncols,nrows,nignore,c,hole
    logical::inside
    
    
    rinc=abs(up-down)/dble(nrows)
    cinc=abs(left-right)/dble(ncols)
    
    nignore=0
    c=0
    do i=1,nrows
        do j=1,ncols
            x=dble(j)*cinc
            y=dble(i)*rinc
            inside=.false.
            do hole=1,nholes
                if(((x-h(hole))**2 + (y-k(hole))**2) < r(hole)) inside=.true.
            end do
            if(.not.inside.and.j<ncols.and.i<nrows)then
                c=c+1
                xynodes(1,c)=x
                xynodes(2,c)=y
            else
                nignore=nignore+1
            end if
        end do
    end do
end
