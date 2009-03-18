subroutine jacobian_fem_jb ( node_num, node_xy, element_num, &
  element_node, quad_num, node_u_variable, node_v_variable, &
  node_p_variable, variable_num, nu, c, ib, a )

!*****************************************************************************80
!
!! JACOBIAN_FEM evaluates the Navier Stokes jacobian matrix.
!
!  Discussion:
!
!    The matrix is known to be banded.  A special matrix storage format
!    is used to reduce the space required.  Details of this format are
!    discussed in the routine DGB_FA.
!
!    The Navier Stokes equations in weak form are:
!
!      Integral ( nu * ( dBdx(I) * dUdx + dBdy(I) * dUdy )
!        + B(I) * ( ( U * dUdx + V * dUdy ) + dPdx - U_RHS ) ) = 0
!
!      Integral ( nu * ( dBdx(I) * dVdx + dBdy(I) * dVdy )
!        + B(I) * ( ( U * dVdx + V * dVdy ) + dPdy - V_RHS ) ) = 0
!
!      Integral ( Q(I) * ( dUdx + dVdy - P_RHS ) ) = 0
!
!    This routine sets up the matrix as though every degree of freedom
!    were unconstrained.  Adjustments for boundary conditions and other
!    constraints should be made after calling this routine.
!
!  Modified:
!
!    08 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer QUAD_NUM, the number of quadrature points used in assembly.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Input, real ( kind = 8 ) C(VARIABLE_NUM), the finite element
!    coefficients of an approximate solution of the Navier Stokes equations.
!
!    Input, integer IB, the bandwidth of the jacobian.
!
!    Output, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the VARIABLE_NUM
!    by VARIABLE_NUM Navier Stokes jacobian, stored in a general band
!    matrix format.
!
  implicit none

  integer ib
  integer node_num
  integer quad_num
  integer element_num
  integer variable_num

  real ( kind = 8 ) a(3*ib+1,variable_num)
  real ( kind = 8 ) area
  real ( kind = 8 ) b(6,quad_num)
  real ( kind = 8 ), dimension(node_num) :: c
  real ( kind = 8 ) cp(3)
  real ( kind = 8 ) cu(6)
  real ( kind = 8 ) cv(6)
  real ( kind = 8 ) dbdx(6,quad_num)
  real ( kind = 8 ) dbdy(6,quad_num)
  real ( kind = 8 ) dpdx(quad_num)
  real ( kind = 8 ) dpdy(quad_num)
  real ( kind = 8 ) dqdx(3,quad_num)
  real ( kind = 8 ) dqdy(3,quad_num)
  real ( kind = 8 ) dudx(quad_num)
  real ( kind = 8 ) dudy(quad_num)
  real ( kind = 8 ) dvdx(quad_num)
  real ( kind = 8 ) dvdy(quad_num)
  integer element
  integer, dimension(6,element_num) :: element_node
  integer i
  integer ip(3)
  integer iu(6)
  integer iv(6)
  integer j
  integer node_p_variable(node_num)
  integer node_u_variable(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) nu
  real ( kind = 8 ) p(quad_num)
  real ( kind = 8 ) p_rhs(quad_num)
  real ( kind = 8 ) q(3,quad_num)
  integer quad
  real ( kind = 8 ), dimension(quad_num) :: quad_w
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy
  real ( kind = 8 ), dimension(2,3) :: t3
  real ( kind = 8 ), dimension(2,6) :: t6
  integer test
  real ( kind = 8 ) triangle_area_2d
  real ( kind = 8 ) u(quad_num)
  real ( kind = 8 ) u_rhs(quad_num)
  real ( kind = 8 ) v(quad_num)
  real ( kind = 8 ) v_rhs(quad_num)
  real ( kind = 8 ) w(quad_num)
  real ( kind = 8 ), dimension(2,quad_num) :: xy
!
!  Initialize the jacobian to zero.
!
  a(1:3*ib+1,1:variable_num) = 0.0D+00
!
!  Get the quadrature weights and nodes.
!
  call quad_rule ( quad_num, quad_w, quad_xy )
!
!  Consider all quantities associated with a given ELEMENT.
!
  do element = 1, element_num
!
!  Extract the nodes of the linear and quadratic triangles.
!
    t3(1:2,1:3) = node_xy(1:2,element_node(1:3,element))
    t6(1:2,1:6) = node_xy(1:2,element_node(1:6,element))
!
!  Map the quadrature points QUAD_XY to points XY in the physical element.
!
    call reference_to_physical_t6 ( t6, quad_num, quad_xy, xy )
    area = abs ( triangle_area_2d ( t3 ) )
    w(1:quad_num) = quad_w(1:quad_num) * area
!
!  Evaluate the basis functions at the quadrature points.
!
    call basis_mn_t6 ( t6, quad_num, xy, b, dbdx, dbdy )

    call basis_mn_t3 ( t3, quad_num, xy, q, dqdx, dqdy )
!
!  Extract the indices of the finite element coefficients for this element.
!
    iu(1:6) = node_u_variable(element_node(1:6,element))
    iv(1:6) = node_v_variable(element_node(1:6,element))
    ip(1:3) = node_p_variable(element_node(1:3,element))
!
!  Extract the finite element coefficients for this element.
!
    cu(1:6) = c(iu(1:6))
    cv(1:6) = c(iv(1:6))
    cp(1:3) = c(ip(1:3))
!
!  Evaluate the flowfield at each quadrature point.
!
     u(1:quad_num)   = matmul ( cu(1:6),  b(1:6,1:quad_num) )
    dudx(1:quad_num) = matmul ( cu(1:6), dbdx(1:6,1:quad_num) )
    dudy(1:quad_num) = matmul ( cu(1:6), dbdy(1:6,1:quad_num) )

     v(1:quad_num)   = matmul ( cv(1:6),  b(1:6,1:quad_num) )
    dvdx(1:quad_num) = matmul ( cv(1:6), dbdx(1:6,1:quad_num) )
    dvdy(1:quad_num) = matmul ( cv(1:6), dbdy(1:6,1:quad_num) )

     p(1:quad_num)   = matmul ( cp(1:3),  q(1:3,1:quad_num) )
    dpdx(1:quad_num) = matmul ( cp(1:3), dqdx(1:3,1:quad_num) )
    dpdy(1:quad_num) = matmul ( cp(1:3), dqdy(1:3,1:quad_num) )
!
!  dUeqn/dUcof,
!  dUeqn/dVcof,
!  dUeqn/dPcof.
!
    do i = 1, 6
      do j = 1, 6

        a(iu(i)-iu(j)+2*ib+1,iu(j)) = a(iu(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) *                                               &
          (                                                             &
            nu * ( dbdx(j,1:quad_num) * dbdx(i,1:quad_num)              &
                 + dbdy(j,1:quad_num) * dbdy(i,1:quad_num) )            &
            +                                                           &
            ( b(j,1:quad_num) * dudx(1:quad_num)                        &
            + u(1:quad_num) * dbdx(j,1:quad_num)                        &
            + v(1:quad_num) * dbdy(j,1:quad_num) ) * b(i,1:quad_num)    &
          )                                                             &
        )

        a(iu(i)-iv(j)+2*ib+1,iv(j)) = a(iu(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) * b(j,1:quad_num) * dudy(1:quad_num)            &
        * b(i,1:quad_num) )

      end do

      do j = 1, 3
        a(iu(i)-ip(j)+2*ib+1,ip(j)) = a(iu(i)-ip(j)+2*ib+1,ip(j)) + sum &
          ( w(1:quad_num) * dqdx(j,1:quad_num) * b(i,1:quad_num) )
      end do

    end do
!
!  dVeqn/dUcof,
!  dVeqn/dVcof,
!  dVeqn/dPcof.
!
    do i = 1, 6
      do j = 1, 6

        a(iv(i)-iu(j)+2*ib+1,iu(j)) = a(iv(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) * b(j,1:quad_num) * dvdx(1:quad_num)            &
          * b(i,1:quad_num) )

        a(iv(i)-iv(j)+2*ib+1,iv(j)) = a(iv(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) *                                               &
          (                                                             &
            nu * ( dbdx(j,1:quad_num) * dbdx(i,1:quad_num)              &
                 + dbdy(j,1:quad_num) * dbdy(i,1:quad_num) )            &
            +                                                           &
            ( u(1:quad_num) * dbdx(j,1:quad_num)                        &
            + b(j,1:quad_num) * dvdy(1:quad_num)                        &
            + v(1:quad_num) * dbdy(j,1:quad_num) )                      &
            * b(i,1:quad_num)                                           &
          )                                                             &
        )

      end do

      do j = 1, 3
        a(iv(i)-ip(j)+2*ib+1,ip(j)) = a(iv(i)-ip(j)+2*ib+1,ip(j)) + sum &
        ( w(1:quad_num) * dqdy(j,1:quad_num) * b(i,1:quad_num) )
      end do

    end do
!
!  dPeqn/dUcof,
!  dPeqn/dVcof,
!
    do i = 1, 3
      do j = 1, 6

        a(ip(i)-iu(j)+2*ib+1,iu(j)) = a(ip(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) * dbdx(j,1:quad_num) * q(i,1:quad_num) )

        a(ip(i)-iv(j)+2*ib+1,iv(j)) = a(ip(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) * dbdy(j,1:quad_num) * q(i,1:quad_num) )

      end do
    end do

  end do

  return
end
