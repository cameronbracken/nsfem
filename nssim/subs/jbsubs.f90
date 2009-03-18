subroutine basis_mn_t3 ( t, n, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_MN_T3: all bases at N points for a T3 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices of a triangle.
!    It works directly with these coordinates, and does not refer to a
!    reference element.
!
!    The sides of the triangle DO NOT have to lie along a coordinate
!    axis.
!
!    The routine evaluates the basis functions associated with each vertex,
!    and their derivatives with respect to X and Y.
!
!  Physical Element T3:
!
!            3
!           / \
!          /   \
!         /     \
!        /       \
!       1---------2
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the vertices
!    of the triangle.  It is common to list these points in counter clockwise
!    order.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the points where the basis functions
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(3,N), the value of the basis functions
!    at the evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(3,N), DPHIDY(3,N), the value of the
!    derivatives at the evaluation points.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) AREA, is (twice) the area of the triangle.
!
  implicit none

  integer n

  real ( kind = 8 ) area
  real ( kind = 8 ) dphidx(3,n)
  real ( kind = 8 ) dphidy(3,n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ) phi(3,n)
  real ( kind = 8 ) t(2,3)

  area = t(1,1) * ( t(2,2) - t(2,3) ) &
       + t(1,2) * ( t(2,3) - t(2,1) ) &
       + t(1,3) * ( t(2,1) - t(2,2) )

  if ( area == 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_MN_T3 - Fatal error!'
    write ( *, '(a)' ) '  Element has zero area.'
    stop

  end if

  phi(1,1:n) =     (   ( t(1,3) - t(1,2) ) * ( p(2,1:n) - t(2,2) )     &
                     - ( t(2,3) - t(2,2) ) * ( p(1,1:n) - t(1,2) ) )
  dphidx(1,1:n) =    - ( t(2,3) - t(2,2) )
  dphidy(1,1:n) =      ( t(1,3) - t(1,2) )

  phi(2,1:n) =     (   ( t(1,1) - t(1,3) ) * ( p(2,1:n) - t(2,3) )     &
                     - ( t(2,1) - t(2,3) ) * ( p(1,1:n) - t(1,3) ) )
  dphidx(2,1:n) =    - ( t(2,1) - t(2,3) )
  dphidy(2,1:n) =      ( t(1,1) - t(1,3) )

  phi(3,1:n) =     (   ( t(1,2) - t(1,1) ) * ( p(2,1:n) - t(2,1) )     &
                     - ( t(2,2) - t(2,1) ) * ( p(1,1:n) - t(1,1) ) )
  dphidx(3,1:n) =    - ( t(2,2) - t(2,1) )
  dphidy(3,1:n) =      ( t(1,2) - t(1,1) )
!
!  Normalize.
!
  phi(1:3,1:n) = phi(1:3,1:n) / area
  dphidx(1:3,1:n) = dphidx(1:3,1:n) / area
  dphidy(1:3,1:n) = dphidy(1:3,1:n) / area

  return
end
subroutine basis_mn_t6 ( t, n, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_MN_T6: all bases at N points for a T6 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices and midside
!    nodes of a triangle.  It works directly with these coordinates, and does
!    not refer to a reference element.
!
!    This routine requires that the midside nodes be "in line"
!    with the vertices, that is, that the sides of the triangle be
!    straight.  However, the midside nodes do not actually have to
!    be halfway along the side of the triangle.
!
!  Physical element T6:
!
!    This picture indicates the assumed ordering of the six nodes
!    of the triangle.
!
!             3
!            / \
!           /   \
!          6     5
!         /       \
!        /         \
!       1-----4-----2
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the nodal oordinates of the element.
!    It is common to list these points in counter clockwise order.
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the coordinates of the point where
!    the basis functions are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(6,N), the basis functions at the
!    evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(6,N), DPHIDY(6,N), the derivatives
!    of the basis functions at the evaluation points.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) AREA, is (twice) the area of the triangle.
!
  implicit none

  integer n

  real ( kind = 8 ) dphidx(6,n)
  real ( kind = 8 ) dphidy(6,n)
  real ( kind = 8 ) gn(n)
  real ( kind = 8 ) gx(n)
  real ( kind = 8 ) hn(n)
  real ( kind = 8 ) hx(n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ) phi(6,n)
  real ( kind = 8 ) t(2,6)
!write(*,*)"the triangle points"
!call print_double_mat(t,2,6)
!write(*,*)"the eval point"
!call print_double_mat(p,2,1)
!
!  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  gn(1:n) = ( t(1,1)   - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( t(2,1)   - t(2,2) )

  hx(1:n) = ( p(1,1:n) - t(1,4) ) * ( t(2,6)   - t(2,4) ) &
          - ( t(1,6)   - t(1,4) ) * ( p(2,1:n) - t(2,4) )

  hn(1:n) = ( t(1,1)   - t(1,4) ) * ( t(2,6)   - t(2,4) ) &
          - ( t(1,6)   - t(1,4) ) * ( t(2,1)   - t(2,4) )

  phi(1,1:n) =     ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(1,1:n) =  (      ( t(2,3) - t(2,2) ) * hx(1:n) &
                   + gx(1:n) * ( t(2,6) - t(2,4) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(1,1:n) = -(      ( t(1,3) - t(1,2) ) * hx(1:n) &
                   + gx(1:n) * ( t(1,6) - t(1,4) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( p(2,1:n) - t(2,1) )

  gn(1:n) = ( t(1,2)   - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( t(2,2)   - t(2,1) )

  hx(1:n) = ( p(1,1:n) - t(1,5) ) * ( t(2,4)   - t(2,5) ) &
          - ( t(1,4)   - t(1,5) ) * ( p(2,1:n) - t(2,5) )

  hn(1:n) = ( t(1,2)   - t(1,5) ) * ( t(2,4)   - t(2,5) ) &
          - ( t(1,4)   - t(1,5) ) * ( t(2,2)   - t(2,5) )

  phi(2,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(2,1:n) =  (      ( t(2,3) - t(2,1) ) * hx(1:n) &
               + gx(1:n) * ( t(2,4) - t(2,5) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(2,1:n) = -(      ( t(1,3) - t(1,1) ) * hx(1:n) &
               + gx(1:n) * ( t(1,4) - t(1,5) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  gn(1:n) = ( t(1,3)   - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( t(2,3)   - t(2,2) )

  hx(1:n) = ( p(1,1:n) - t(1,6) ) * ( t(2,5)   - t(2,6) ) &
          - ( t(1,5)   - t(1,6) ) * ( p(2,1:n) - t(2,6) )

  hn(1:n) = ( t(1,3)   - t(1,6) ) * ( t(2,5)   - t(2,6) ) &
          - ( t(1,5)   - t(1,6) ) * ( t(2,3)   - t(2,6) )

  phi(3,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(3,1:n) =  (      ( t(2,1) - t(2,2) ) * hx(1:n) &
               + gx(1:n) * ( t(2,5) - t(2,6) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(3,1:n) = -(      ( t(1,1) - t(1,2) ) * hx(1:n) &
               + gx(1:n) * ( t(1,5) - t(1,6) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,3) ) * ( t(2,1)   - t(2,3) ) &
          - ( t(1,1)   - t(1,3) ) * ( p(2,1:n) - t(2,3) )

  gn(1:n) = ( t(1,4)   - t(1,3) ) * ( t(2,1)   - t(2,3) ) &
          - ( t(1,1)   - t(1,3) ) * ( t(2,4)   - t(2,3) )

  hx(1:n) = ( p(1,1:n) - t(1,3) ) * ( t(2,2)   - t(2,3) ) &
          - ( t(1,2)   - t(1,3) ) * ( p(2,1:n) - t(2,3) )

  hn(1:n) = ( t(1,4)   - t(1,3) ) * ( t(2,2)   - t(2,3) ) &
          - ( t(1,2)   - t(1,3) ) * ( t(2,4)   - t(2,3) )

  phi(4,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(4,1:n) =  (      ( t(2,1) - t(2,3) ) * hx(1:n) &
               + gx(1:n) * ( t(2,2) - t(2,3) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(4,1:n) = -(      ( t(1,1) - t(1,3) ) * hx(1:n) &
               + gx(1:n) * ( t(1,2) - t(1,3) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,1) ) * ( t(2,2)   - t(2,1) ) &
          - ( t(1,2)   - t(1,1) ) * ( p(2,1:n) - t(2,1) )

  gn(1:n) = ( t(1,5)   - t(1,1) ) * ( t(2,2)   - t(2,1) ) &
          - ( t(1,2)   - t(1,1) ) * ( t(2,5)   - t(2,1) )

  hx(1:n) = ( p(1,1:n) - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( p(2,1:n) - t(2,1) )

  hn(1:n) = ( t(1,5)   - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( t(2,5)   - t(2,1) )

  phi(5,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(5,1:n) =  (      ( t(2,2) - t(2,1) ) * hx(1:n) &
               + gx(1:n) * ( t(2,3) - t(2,1) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(5,1:n) = -(      ( t(1,2) - t(1,1) ) * hx(1:n) &
               + gx(1:n) * ( t(1,3) - t(1,1) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  gn(1:n) = ( t(1,6)   - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( t(2,6)   - t(2,2) )

  hx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  hn(1:n) = ( t(1,6)   - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( t(2,6)   - t(2,2) )

  phi(6,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(6,1:n) =  (      ( t(2,1) - t(2,2) ) * hx(1:n) &
               + gx(1:n) * ( t(2,3) - t(2,2) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(6,1:n) = -(      ( t(1,1) - t(1,2) ) * hx(1:n) &
               + gx(1:n) * ( t(1,3) - t(1,2) ) ) / ( gn(1:n) * hn(1:n) )
  	!write(*,*)'one of the basis functions'
	!call print_double_mat(phi,6,1)
  return
end
subroutine quad_rule ( quad_num, quad_w, quad_xy )

!*****************************************************************************80
!
!! QUAD_RULE sets the quadrature rule for assembly.
!
!  Discussion:
!
!    The quadrature rule is given for a reference element.
!
!      0 <= X,
!      0 <= Y, and
!      X + Y <= 1.
!
!      ^
!    1 | *
!      | |\
!    Y | | \
!      | |  \
!    0 | *---*
!      +------->
!        0 X 1
!
!    The rules have the following precision:
!
!    QUAD_NUM  Precision
!
!     1        1
!     3        2
!     4        3
!     6        4
!     7        5
!     9        6
!    13        7
!
!  Modified:
!
!    18 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer QUAD_NUM, the number of quadrature nodes.
!
!    Output, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Output, real ( kind = 8 ) QUAD_XY(2,QUAD_NUM),
!    the coordinates of the quadrature nodes.
!
  implicit none

  integer quad_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ), dimension(quad_num) :: quad_w
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w

  if ( quad_num == 1 ) then

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      1.0D+00 / 3.0D+00, 1.0D+00 / 3.0D+00 /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = 1.0D+00

  else if ( quad_num == 3 ) then

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      0.5D+00, 0.0D+00, &
      0.5D+00, 0.5D+00, &
      0.0D+00, 0.5D+00 /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = 1.0D+00 / 3.0D+00

  else if ( quad_num == 4 ) then

    a =   6.0D+00 / 30.0D+00
    b =  10.0D+00 / 30.0D+00
    c =  18.0D+00 / 30.0D+00

    d =  25.0D+00 / 48.0D+00
    e = -27.0D+00 / 48.0D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      b, b, &
      c, a, &
      a, c, &
      a, a /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ e, d, d, d /)

  else if ( quad_num == 6 ) then

    a = 0.816847572980459D+00
    b = 0.091576213509771D+00
    c = 0.108103018168070D+00
    d = 0.445948490915965D+00
    v = 0.109951743655322D+00
    w = 0.223381589678011D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      a, b, &
      b, a, &
      b, b, &
      c, d, &
      d, c, &
      d, d /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ v, v, v, w, w, w /)

  else if ( quad_num == 7 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -       sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +       sqrt ( 15.0D+00 ) ) / 21.0D+00
    u = 0.225D+00
    v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
    w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      a, a, &
      b, c, &
      c, b, &
      c, c, &
      d, e, &
      e, d, &
      e, e /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ u, v, v, v, w, w, w /)

  else if ( quad_num == 9 ) then

    a = 0.124949503233232D+00
    b = 0.437525248383384D+00
    c = 0.797112651860071D+00
    d = 0.165409927389841D+00
    e = 0.037477420750088D+00

    u = 0.205950504760887D+00
    v = 0.063691414286223D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      a, b, &
      b, a, &
      b, b, &
      c, d, &
      c, e, &
      d, c, &
      d, e, &
      e, c, &
      e, d /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ u, u, u, v, v, v, v, v, v /)

  else if ( quad_num == 13 ) then

    h = 1.0D+00 / 3.0D+00
    a = 0.479308067841923D+00
    b = 0.260345966079038D+00
    c = 0.869739794195568D+00
    d = 0.065130102902216D+00
    e = 0.638444188569809D+00
    f = 0.312865496004875D+00
    g = 0.048690315425316D+00

    w = -0.149570044467670D+00
    t =  0.175615257433204D+00
    u =  0.053347235608839D+00
    v =  0.077113760890257D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      h, h, &
      a, b, &
      b, a, &
      b, b, &
      c, d, &
      d, c, &
      d, d, &
      e, f, &
      e, g, &
      f, e, &
      f, g, &
      g, e, &
      g, f /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ w, t, t, t, u, u, u, v, v, v, v, v, v /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUAD_RULE - Fatal error!'
    write ( *, '(a,i8)' ) '  No rule is available of order QUAD_NUM = ', &
      quad_num
    stop

  end if

  return
end
subroutine reference_to_physical_t6 ( t, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 6 physical triangle and a point
!    (XSI,ETA) in the reference triangle, the routine computes the value
!    of the corresponding image point (X,Y) in physical space.
!
!    The mapping from (XSI,ETA) to (X,Y) has the form:
!
!      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!                 + D1 * XSI    + E1 * ETA     + F1
!
!      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!                 + D2 * XSI    + E2 * ETA     + F2
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Modified:
!
!    25 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0), (1,0),
!    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
!
!    Input, integer N, the number of objects to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference triangle.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical triangle.
!
  implicit none

  integer n

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  real ( kind = 8 ) c(2)
  real ( kind = 8 ) d(2)
  real ( kind = 8 ) e(2)
  real ( kind = 8 ) f(2)
  integer i
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,6)

  do i = 1, 2

    a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2)                    &
           - 4.0D+00 * t(i,4)

    b(i) =   4.0D+00 * t(i,1)                                       &
           - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

    c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3) &
                                                 - 4.0D+00 * t(i,6)

    d(i) = - 3.0D+00 * t(i,1) -           t(i,2)                    &
           + 4.0D+00 * t(i,4)

    e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3) &
                                                 + 4.0D+00 * t(i,6)
    f(i) =             t(i,1)

  end do

  do i = 1, 2
    phy(i,1:n) = a(i) * ref(1,1:n) * ref(1,1:n) &
               + b(i) * ref(1,1:n) * ref(2,1:n) &
               + c(i) * ref(2,1:n) * ref(2,1:n) &
               + d(i) * ref(1,1:n) &
               + e(i) * ref(2,1:n) &
               + f(i)
  end do

  return
end
function triangle_area_2d ( t )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) TRIANGLE_AREA_2D, the absolute area of
!    the triangle.
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) triangle_area_2d

  triangle_area_2d = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end


subroutine bandwidth ( element_order, element_num, element_node, &
  node_num, node_p_variable, node_u_variable, node_v_variable, ib )

!*****************************************************************************80
!
!! BANDWIDTH determines the bandwidth of the coefficient matrix.
!
!  Discussion:
!
!    We take the bandwidth to be the maximum difference between the
!    indices of two variables associated with nodes that share an element.
!
!    Therefore, we can compute the bandwidth by examining each element,
!    and finding the maximum difference in indices of any two variables
!    associated with nodes in that element.
!
!  Modified:
!
!    10 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_ORDER, the number of nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Output, integer IB, the half bandwidth of the matrix.
!
  implicit none

  integer node_num
  integer element_num
  integer element_order

  integer element
  integer element_node(element_order,element_num)
  integer i4_huge
  integer local
  integer ib
  integer node
  integer node_p_variable(node_num)
  integer node_u_variable(node_num)
  integer node_v_variable(node_num)
  integer v
  integer v_max
  integer v_min

  ib = 0

  do element = 1, element_num

    v_max = -i4_huge ( )
    v_min = i4_huge ( )

    do local = 1, element_order

      node = element_node(local,element)

      v = node_u_variable(node)
      v_max = max ( v_max, v )
      v_min = min ( v_min, v )

      v = node_v_variable(node)
      v_max = max ( v_max, v )
      v_min = min ( v_min, v )

      if ( 0 < node_p_variable(node) ) then
        v = node_p_variable(node)
        v_max = max ( v_max, v )
        v_min = min ( v_min, v )
      end if

    end do

    ib = max ( ib, v_max - v_min )

  end do

  return
end

subroutine dgb_fa ( n, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! DGB_FA performs a LINPACK-style PLU factorization of an DGB matrix.
!
!  Discussion:
!
!    The DGB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!    The following program segment will set up the input.
!
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i - j + m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of the array A.
!    In addition, the first ML rows in the array are used for
!    elements generated during the triangularization.
!
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Modified:
!
!    04 March 1999
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User''s Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real ( kind = 8 ) A(2*ML+MU+1,N), on input, the matrix
!    in band storage, on output, information about the LU factorization.
!
!    Output, integer PIVOT(N), the pivot vector.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ml
  integer mu
  integer n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer i
  integer i0
  integer info
  integer pivot(n)
  integer j
  integer j0
  integer j1
  integer ju
  integer jz
  integer k
  integer l
  integer lm
  integer m
  integer mm
  real ( kind = 8 ) t
  real ( kind = 8 ) temp

  m = ml + mu + 1
  info = 0
!
!  Zero out the initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    a(i0:ml,jz) = 0.0D+00
  end do

  jz = j1
  ju = 0

  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      a(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )

    l = m
    do j = m+1, m+lm
      if ( abs ( a(l,k) ) < abs ( a(j,k) ) ) then
        l = j
      end if
    end do

    pivot(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange if necessary.
!
    temp = a(l,k)
    a(l,k) = a(m,k)
    a(m,k) = temp
!
!  Compute multipliers.
!
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  Row elimination with column indexing.
!
    ju = max ( ju, mu+pivot(k) )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju

      l = l - 1
      mm = mm - 1

      if ( l /= mm ) then
        temp = a(l,j)
        a(l,j) = a(mm,j)
        a(mm,j) = temp
      end if

      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)

    end do

  end do

  pivot(n) = n
  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGB_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end

subroutine dgb_sl ( n, ml, mu, a, pivot, b, job )

!*****************************************************************************80
!
!! DGB_SL solves a system factored by DGB_FA.
!
!  Discussion:
!
!    The DGB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!  Modified:
!
!    04 March 1999
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User''s Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the LU factors from DGB_FA.
!
!    Input, integer PIVOT(N), the pivot vector from DGB_FA.
!
!    Input/output, real ( kind = 8 )l B(N).
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Input, integer JOB.
!    0, solve A * x = b.
!    nonzero, solve A'' * x = b.
!
  implicit none

  integer ml
  integer mu
  integer n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer pivot(n)
  integer j
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    if ( 1 <= ml ) then

      do k = 1, n-1

        lm = min ( ml, n-k )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a(m+1:m+lm,k)

      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a(la:la+lm-1,k)

    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      b(k) = ( b(k) - sum ( a(la:la+lm-1,k) * b(lb:lb+lm-1) ) ) &
        / a(m,k)
    end do
!
!  Solve L'' * X = Y.
!
    if ( 1 <= ml ) then

      do k = n-1, 1, -1

        lm = min ( ml, n-k )
        b(k) = b(k) + sum ( a(m+1:m+lm,k) * b(k+1:k+lm) )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

      end do

    end if

  end if

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer I4_HUGE, a "huge" integer.
!
  implicit none

  integer i4_huge

  i4_huge = huge ( 1 )

  return
end

subroutine jacobian_adjust_dirichlet ( node_num, node_xy, onboundary,node_u_variable, &
  node_v_variable, node_p_variable, variable_num, ib, a )

!*****************************************************************************80
!
!! JACOBIAN_ADJUST_DIRICHLET adjusts the jacobian for Dirichlet conditions.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, integer IB, the half-bandwidth of the matrix.
!
!    Input/output, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the VARIABLE_NUM
!    by VARIABLE_NUM coefficient matrix, stored in a compressed format;
!    on output, the matrix has been adjusted for Dirichlet boundary conditions.
!
  implicit none

  integer ib
  integer node_num
  integer variable_num

  real ( kind = 8 ), dimension(3*ib+1,variable_num) :: a
  integer column
  integer column_high
  integer column_low
  integer, parameter :: DIRICHLET = 2
  integer ip
  integer iu
  integer iv
  integer node
  integer onboundary(variable_num)
  integer node_p_variable(node_num)
  !  integer node_u_condition(node_num)
  integer node_u_variable(node_num)
  !  integer node_v_condition(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ) node_xy(2,node_num)
  logical::pressure
     pressure=.false.
 do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( onboundary(node) >= 1 ) then

      column_low = max ( iu - ib, 1 )
      column_high = min ( iu + ib, variable_num )

      do column = column_low, column_high
        a(iu-column+2*ib+1,column) = 0.0D+00
      end do
      a(2*ib+1,iu) = 1.0D+00

    end if

    if ( onboundary(node) >= 1 ) then

      column_low = max ( iv - ib, 1 )
      column_high = min ( iv + ib, variable_num )

      do column = column_low, column_high
        a(iv-column+2*ib+1,column) = 0.0D+00
      end do
      a(2*ib+1,iv) = 1.0D+00

    end if

    if ( 0 < ip ) then
      if ( onboundary(node)==1.and..not.pressure.and.node_p_variable(node) >= 0 ) then
        pressure=.true.
        write(*,*)'Specifying pressure at node',node
        column_low = max ( ip - ib, 1 )
        column_high = min ( ip + ib, variable_num )

        do column = column_low, column_high
          a(ip-column+2*ib+1,column) = 0.0D+00
        end do
        a(2*ib+1,ip) = 1.0D+00

      end if

    end if

  end do

  return
end
