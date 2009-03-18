subroutine boundary_conditions(A,onboundary,n,nvars,u_position,v_position,&
                                                        p_position,filename)
	implicit none
	integer::n,i,nvars
    character(len=*)::filename
	double precision,dimension(nvars,nvars),intent(out)::A
	integer,dimension(n)::onboundary,u_position,v_position,p_position
	logical::pressure,verbose

	!!This subroutine imposes boundary conditions on the coefficient matrix A
	!!It is the counterpart to the boundary_type subtroutine but it is general 
	!!and does not need to be modified for a specific probelm.
	!!The routine imposes dirchlet conditions on the coefficient matrix A by zeroing
	!!out a row associated with a boundary and setting A(i,i)=1 effectively 
	!!forcing the noe to have a specific value.
	!!
	!!NOTE: Neumann conditions are not implemented, though they are easier to handle.
	!!
	!!Variables(inputs):
    !!n             - [integer] number of nodes
	!!nvars         - [integer] number of variables or equations (not n*3 because linear
    !!                           pressure is not defined at midside nodes)
	!!A				- [double (nvars,nvars)] the coefficient matrix, which is modified 
	!!                                       and passed back 
    !!node_xy       - [double (n,2)] coordinates of each node
    !!onboundary    - [integer (n)] 1 if on boundary, 0 if not 
    !!f             - [double (nvars)] the right hand side vector
    !!u_position    - [integer (n)] the indicies of the horizontal velocity equations
    !!v_position    - [integer (n)] the indicies of the vertical velocity equations
    !!p_position    - [integer (n)] the indicies of the pressure equations, if -1 then 
    !!                              the node is a midside node and pressure is not defined


	open(15,file=filename)
	pressure=.false.
	verbose=.false.

	do i=1,n

		if(onboundary(i)>0)then 
			A(u_position(i),:)=0
			A(v_position(i),:)=0
			A(u_position(i),u_position(i))=1d0
			A(v_position(i),v_position(i))=1d0
		end if
    
		if(.not.pressure.and.p_position(i)>0)then
            !write(*,'(a,i3,x,a,x,2f10.3)')'Specifying pressure variable',p_position(i)
		    A(p_position(i),:)=0
            A(p_position(i),p_position(i))=1
			pressure=.true.
		end if

		write(15,*)onboundary(i)

	end do

	close(15)
	return
end
