subroutine boundary_type(node_xy,onboundary,f,n,nvars,u_position,v_position,p_position)
    implicit none
    integer::i,n,nvars
    integer,dimension(n),intent(inout)::onboundary
    integer,dimension(n),intent(in)::u_position,v_position,p_position
    double precision,dimension(n,2),intent(in)::node_xy
    double precision::up,down,left,right,h,k,ninteriornodes,nsidepts,nobspts,radius
    double precision::in,th,xi,yi,tol
    logical::verbose,pressure,hole
    character(len=10)::oshape
    double precision,dimension(nvars)::f
    
    !!This subroutine determines the type of the boundary condition at every node,
    !!It is intended to be problem specific. The vector onboundary is passed in 
    !!already containing which nodes are on ANY boundary of the system but this 
    !!routine is intended to classify the interior boundary nodes and pass onboundary
    !!back out:
    !!  0=not on boundary
    !!  1=on exterior boundary
    !!  2=on an interior boundary (hole)
    !!
    !!This subroutine also looks for a file called 'domain' which for this problem 
    !!specifies the dimensions of the problem domain and the hole
    !!
    !!
    !!Variables(inputs):
    !!n             - [integer] number of nodes
    !!nvars         - [integer] number of variables or equations (not n*3 because linear
    !!                          pressure is not defined at midside nodes)
    !!node_xy       - [double (n,2)] coordinates of each node
    !!onboundary    - [integer (n)] 1 if on boundary, 0 if not 
    !!f             - [double (nvars)] the right hand side vector
    !!u_position    - [integer (n)] the indicies of the horizontal velocity equations
    !!v_position    - [integer (n)] the indicies of the vertical velocity equations
    !!p_position    - [integer (n)] the indicies of the pressure equations, if -1 then 
    !!                              the node is a midside node and pressure is not defined
    !!(inputs form the domain file)
    !!up,down,left,right - [double] the boundaries of the system
    !!h,k                - [double] midpoint of the obstrction
    !!radius             - [double] radius of the hole
    !!oshape             - [character] the shape of the obstruction
    
        
    open(15,file='domain')
    read(15,*)up,down,left,right
    close(15)
    
    verbose=.false.

    in=2d0

    pressure=.false.
    f=0
    tol=0.0001

    do i=1,n
        xi=node_xy(i,1)
        yi=node_xy(i,2)
        hole=.true.

        if(onboundary(i)==1)then
            if(dabs(xi-left)<tol)then  !on left boundary
                    onboundary(i)=1
                    hole=.false.
            end if 
            if(dabs(xi-right)<tol)then  !on right boundary
                    onboundary(i)=1
                    hole=.false.
            end if
            if(dabs(yi-up)<tol)then  !on top boundary
                    hole=.false.
                if(xi > (.01*dabs(left-right)+left) .and. xi < (.99*dabs(left-right)+left))then
                    onboundary(i)=0
                else
                    onboundary(i)=1
                end if
            end if
            
            if(dabs(yi-down)<tol)then  !on bottom boundary
                onboundary(i)=1
                    hole=.false.
                    f(v_position(i)) = in*xi*(right-xi)/dabs(right-left)
            end if

            if(hole)then   !on interior boundary
                    onboundary(i)=2
            end if
        end if

        if(.not.pressure.and.p_position(i)>0)then
                write(*,'(a,i3,x,a,x,i5,x,a,x,2f10.3)')'Specifying pressure variable',p_position(i),'node',n,'at',node_xy(i,1),node_xy(i,2)
                f(p_position(i))=0
                pressure=.true.
        end if

    end do
    return
end

