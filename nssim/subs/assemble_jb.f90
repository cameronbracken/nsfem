subroutine assemble_jb(A,f,node_xy,element_node,node_num,element_num,&
			quad_num,variable_num,nu,node_u_variable,&
            node_v_variable,node_p_variable)
	implicit none
	integer::node_num,element_num,quad_num,variable_num
	double precision,dimension(2,node_num)::node_xy
	integer,dimension(6,element_num)::element_node
	integer,dimension(node_num)::node_u_variable,node_v_variable,node_p_variable
	double precision::nu
	double precision,dimension(variable_num,variable_num)::A
	double precision,dimension(variable_num)::f
	double precision,dimension(quad_num)::quad_w
	double precision,dimension(2,quad_num)::quad_xy
	double precision,dimension(2,quad_num)::xy
	double precision::triangle_area_2d,area
	double precision,dimension(quad_num)::w
	double precision,dimension(2,1)::point
	double precision,dimension(2,6)::t6
	double precision,dimension(2,3)::t3
	double precision,dimension(6,1)::bi, dbidx, dbidy ,bj, dbjdx, dbjdy 
	double precision,dimension(3,1)::qi, dqidx, dqidy ,qj, dqjdx, dqjdy
	integer::basis_node,test_node
	integer::element,quad,test,basis
	integer::iu,ju,iv,jv,ip,jp


    	f = 0d0
		A = 0d0

		!  Get the quadrature weights and nodes.
  call quad_rule (quad_num,quad_w,quad_xy)


  do element=1,element_num

		!  Extract the nodes of the linear and quadratic triangles.
    t3(1:2,1:3) = node_xy(1:2,element_node(1:3,element))
    t6(1:2,1:6) = node_xy(1:2,element_node(1:6,element))
	
	
	!  Map the quadrature points QUAD_XY to points XY in the physical triangle.
    call reference_to_physical_t6(t6,quad_num,quad_xy,xy)
    area = abs(triangle_area_2d(t3))
    w(1:quad_num) = area*quad_w(1:quad_num);
	![ u_rhs, v_rhs, p_rhs ] = rhs ( quad_num, xy );
	!
	!  Consider the QUAD-th quadrature point.
	!


    do quad = 1,quad_num

      point(:,1) = xy(1:2,quad)
	
	
		!Evaluate the test functions.
      call basis_mn_t6(t6,1,point,bi,dbidx,dbidy)
      call basis_mn_t3(t3,1,point,qi,dqidx,dqidy)

      do test = 1 , 6

        test_node = element_node(test,element);

        iu = node_u_variable(test_node);
        iv = node_v_variable(test_node);
        ip = node_p_variable(test_node);
		
		
		!
		!  Compute the source terms do the right hand side.
		!
        !f(iu) = f(iu) + w(quad) * u_rhs(quad) * bi(test);
        !f(iv) = f(iv) + w(quad) * v_rhs(quad) * bi(test);
        !if ( 0 < ip )
        !  f(ip) = f(ip) + w(quad) * p_rhs(quad) * qi(test);
        !end
		
		
		
		!  Consider the basis functions.
        call basis_mn_t6(t6,1,point,bj,dbjdx,dbjdy)
		call basis_mn_t3(t3,1,point,qj,dqjdx,dqjdy)

        do basis = 1,6

          basis_node = element_node(basis,element);

          ju = node_u_variable(basis_node);
          jv = node_v_variable(basis_node);
          jp = node_p_variable(basis_node);
	

			!  Add terms to the horizonal momentum equation.
          a(iu,ju) = a(iu,ju) + w(quad)*nu*(dbidx(test,1)*dbjdx(basis,1)+dbidy(test,1)*dbjdy(basis,1))

          if (jp>0)then !we are at a corner node
            a(iu,jp) = a(iu,jp) + w(quad)*bi(test,1)*dqjdx(basis,1)
          end if
		
		
			!  Add terms to the vertical momentum equation.
          a(iv,jv) = a(iv,jv) + w(quad)*nu*(dbidx(test,1)*dbjdx(basis,1)+dbidy(test,1)*dbjdy(basis,1))

          if(jp>0)then !we are at a corner node
            a(iv,jp) = a(iv,jp) + w(quad)*bi(test,1)*dqjdy(basis,1)
          end if
		
		
			!  Add terms to the continuity equation.
          if(ip>0)then !we are at a corner node
            a(ip,ju) = a(ip,ju) + w(quad)*qi(test,1)*dbjdx(basis,1)
            a(ip,jv) = a(ip,jv) + w(quad)*qi(test,1)*dbjdy(basis,1)
          end if

        end do

      end do

    end do

  end do

return
end subroutine
