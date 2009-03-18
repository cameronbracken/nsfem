program nssim
    implicit none
    character(len=50)::node_file,ele_file 
    double precision,dimension(:,:),allocatable::node_xy,node_xy_t,ele_num_d
    double precision,dimension(:,:),allocatable::A
    double precision,dimension(:),allocatable::f,x,delx,r
    integer,dimension(:,:),allocatable::ele_num,ele_num_t
    integer,dimension(:,:),allocatable::hold
    integer::i,j,nnode,nele,c1,c2,nvars
    logical::verbose,quiet,sparse,ns,jb
    integer,dimension(:),allocatable::onboundary,u_position,v_position,p_position
    integer,dimension(:),allocatable::ipiv
    integer::lda,info,nrhs,ldb
    integer,dimension(:),allocatable::corner_node
    integer::nqp,ncnode,nelt,ns_it,ib
    double precision::nu,nstol
    integer,allocatable,dimension(:)::ia,ja,ia_guess,ja_guess
    double precision,allocatable,dimension(:)::aa,aa_guess
    double precision::OMP_GET_WTIME,time,solvetime,assembletime,alpha

    !!This is the main navier stokes simulation program to solve the steady
    !!state navier-stokes equations.  The program looks for two files called 
    !!data_elements.in and data_nodes.in containing the  element indicies and 
    !!the node coordinates.  The main components of the program involve calling
    !!a routine to assemble the finite element coefficient matrix and finally to
    !!call the LAPACK routines DGETRF and DGETRS to solve the system of linear
    !!equations.
    !!
    !!~Input Files~
    !!data_elments.in - is intended to be the output .ele file from the 
    !!                  meshing software triangle, with the first line of the
    !!                  default .ele file removed as if created with the command
    !!                      sed '1d' mesh.1.ele > data_elements.in
    !!                  as such it will have 7 colums, the first one being the
    !!                  element number followed by the indicies of the standard
    !!                  2nd order element (6 nodes per triangle, with the -o2
    !!                  triangle swich). Triangle orders the elements like:
    !!                                  3
    !!                                 / \
    !!                                5   4
    !!                               /     \
    !!                              1---6---2
    !!                  The ordering will be rearranged later in in the program
    !!                  to fit the ordering required by John Burhardt's
    !!                  subroutines.
    !!data_nodes.in   - Is intended to be the output .node file from the meshing
    !!                  software triangle, with the first line of the default
    !!                  .node file removed as if created with the command:
    !!                      sed '1d' mesh.1.node > data_nodes.in
    !!                  The fourth column in file should be the boundary
    !!                  information (see the onboundary variable). 
    !!
    !!~Variables (In order of appearance)~
    !!node_file       - [character] local name of the element input file 
    !!ele_file        - [character] local name of the node input file
    !!verbose         - [logical] should hella crap be printed out
    !!quiet           - [logical] should a reasonable amount of crap be printed
    !!nqp             - [integer] number of gaussian quadriture points
    !!nu              - [double] the kinematic viscocity of the fluid
    !!c1              - [integer] number of columns in the ele file
    !!c2              - [integer] number of columns in the node file
    !!nele            - [integer] number of elements
    !!nnode           - [integer] number of nodes 
    !!ele_num         - [integer (nele,6)] element indicies, rows are elements
    !!node_xy         - [double (nnode,2)] node coordinates
    !!onboundary      - [integer (nnode)] is the node on a boundary 
    !!                          0 - the node is in the interior of the domain
    !!                          1 - the node is on the boundary of the domain 
    !!                          2 - the node is on an interior boundary
    !!                  the .node file from triangle tells us 0 or 1 if the node 
    !!                  is on ANY boundary but the routine boundary_type
    !!                  determines the interior boundary part (for holes)
    !!corner_node     - [integer (nnode)] is the node a corner (vertex) node
    !!                          0 - the node is a midside node 
    !!                          1 - the node is a corner (vertex) node
    !!ncnode          - [integer] the number of corner (vertex) nodes
    !!u_position      - [integer (nnode)] the indicies of the horizontal 
    !!                  velocity equations
    !!v_position      - [integer (nnode)] the indicies of the vertical velocity
    !!                  equations
    !!p_position      - [integer (nnode)] the indicies of the pressure
    !!                  equations, if -1 then the node is a midside node and 
    !!                  pressure is not defined
    !!nvars           - [integer] number of variables (equations) (not n*3     
    !!                  because linear pressure is not defined at midside 
    !!                  nodes)
    !!A               - [double (nvars,nvars)] The coefficient matrix
    !!f               - [double (nvars)] The right hand side vector
    !!
    !!~Output Files~
    !!summary.out     - Gives basic summary statistics regarding computation
    !!                  time and number of variables
    !!elements.out    - Is the vertex element indicies (only corner nodes)
    !!nodes.out       - Is the nodes and nothing else
    !!velocity.out    - Is the u,v (horizontal,vertical) velocity corresponding
    !!                  to the nodes in nodes.out
    !!pressure.out    - is the pressure at the nodes corresponding to nodes.out
    !!--------------------------------------------------------------------------

    call getarg(1,ele_file)
    call getarg(2,node_file)

    !ele_file='elements.in'
    !node_file='nodes.in'
    
    !Parameters
    open(77,file='simfile')
    read(77,*)verbose,quiet,sparse,nqp,nu
    close(77)
    !verbose=.false.
    !quiet=.true.
    !sparse=.true.
    !nqp=3
    !nu=.1
    ns=.false.
    jb=.false.
    !sparse=.true.
    nstol=.2d0
   
    time=OMP_GET_WTIME()
    
    open(22,file='summary.out')
    open(122,file='time.out')

    write(*,*)
    write(*,*)'Navier-Stokes simulation'
    if(.not.quiet)then
        write(*,*)''    
        write(*,*)'Determining number of elements in ',trim(ele_file)
    end if

    call itable_header_read(ele_file,c1,nele)
    
    write(*,*)'There were',nele,'elements'
    write(22,*)'There were',nele,'elements'
    if(.not.quiet)then
        write(*,*)'There are',nele,'elements in ',trim(ele_file),' and',c1,'columns'
        write(*,*)''
        write(*,*)'Determining number of nodes in ',trim(node_file)
    end if 

    call dtable_header_read(node_file,c2,nnode)
    
    write(22,*)'There were',nnode,'nodes'
    if(.not.quiet)then
        write(*,*)'There are',nnode,'nodes in ',trim(node_file),' and',c2,'columns'
        write(*,*)
    end if 

    if(.not.quiet)then
        write(*,*)'Allocating Arrays...'
        write(*,*)''
    end if
    allocate(ele_num_d(nele,c1),ele_num(nele,6),ele_num_t(6,nele),node_xy(nnode,c2),node_xy_t(2,nnode),onboundary(nnode),&
                corner_node(nnode),u_position(nnode),v_position(nnode),p_position(nnode))

    open(11,file=ele_file)
    open(12,file=node_file)

    if(.not.quiet)then
        write(*,*)'Reading Data from ',trim(ele_file)
    end if

    ele_num=0
    do i=1,nele
        read(11,*)ele_num_d(i,:)
    end do
    do i=1,nele
        do j=2,7
            ele_num(i,j-1)=int(ele_num_d(i,j))
        end do
        if(verbose)then
            write(*,*)'read this:',ele_num(i,:)
        end if
    end do
    close(11)
    if(.not.quiet)then
        write(*,*)'Read',c1*nele,'items'
    
        write(*,*)'Reading Data from ',trim(node_file)
    end if

    do i=1,nnode
        read(12,*)node_xy(i,:)
    end do
    if(.not.quiet)then
        write(*,*)'Read',c2*nnode,'items'
        write(*,*)''
    end if
    close(12)

    if(verbose)then
        write(*,*)'Here is the first bit of the input files, how do they look'  
        write(*,*)'Element file:'
        call print_int_mat(ele_num(1:10,:),10,6)
        write(*,*)'Node file:'
        do i=1,10
            write(*,'(1i3,5f10.3,1i2)')int(node_xy(i,1)),node_xy(i,2:3),int(node_xy(i,4))
        end do
        write(*,*)''
    end if
    
    


        !determine which nodes are the corner nodes of an element 
    if(.not.quiet)then
        write(*,*)'Determining Corner Nodes...'
    end if
    corner_node = 0 
    do i=1,nele
        do j=1,3
            corner_node(ele_num(i,j))=1
        end do
    end do

    ncnode=0
    do i=1,nnode
        if(corner_node(i)==1)then
            ncnode=ncnode+1
        end if
    end do
    if(.not.quiet)then
        write(*,*)'There were',ncnode,'vertex nodes.'
    end if

    if(verbose)then 
        write(*,*)'Corner Nodes:'
        call print_int_mat(corner_node,nnode,1)
        write(*,*)
    end if

    open(16,file='corner.out')
    do i=1,nnode
        write(16,*)corner_node(i)
    end do 
    if(.not.quiet)then
        write(*,*)'Wrote corner node info to corner.out'
        write(*,*)
    end if
    close(16)



    if(.not.quiet)then
        write(*,*)'Determining variable positions' 
    end if
    nvars = 0
!$OMP do ordered
    do i = 1 , nnode
!$OMP ordered  
        nvars = nvars + 1
        u_position(i) = nvars

        nvars = nvars + 1
        v_position(i) = nvars

        if ( corner_node(i) == 1 ) then
            nvars = nvars + 1
            p_position(i) = nvars
        else
            p_position(i) = -1
        end if
!$OMP end ordered
    end do
!$OMP end do

    write(*,*)'There were',nvars,'variables.'
    write(22,*)'There were',nvars,'variables.'

        
    allocate(hold(nele,3))
    
    if(.not.quiet)then
        write(*,*)'Rearranging node file to fit John Burkardts format..'
        !John Burkhardt's subroutines assumes a different midside node ordering 
        !than Triangle outputs so the columns need to be swapped around
    end if
    hold=ele_num(:,4:6)
    ele_num(:,4)=hold(:,3)
    ele_num(:,5)=hold(:,1)
    ele_num(:,6)=hold(:,2)
    deallocate(hold)
    if(verbose)then
        write(*,*)'nodes after rearrangeing'
        call print_int_mat(ele_num,nele,6)

        write(*,*)'Coordinates for input:'
        call print_double_mat(node_xy(1:10,2:3),10,2)
    end if
    if(.not.quiet)then
        write(*,*)'Finished rearrangeing.'
        write(*,*)
    end if
    
    open(17,file='elements.out')
    do i=1,nele
        write(17,*)ele_num(i,1:3)
    end do
    if(.not.quiet)then
        write(*,*)'Wrote three corner element file to elements.out'
        write(*,*)
    end if

    open(15,file='nodes.out')
    do i=1,nnode
        write(15,*)node_xy(i,2),node_xy(i,3)
    end do
    if(.not.quiet)then
        write(*,*)'Wrote node information to nodes.out'
        write(*,*)
    end if
    
    
     
    onboundary(:)=int(node_xy(:,c2))
    
    if(.not.quiet)then
        write(*,*)'Calling assembly subroutine...'
        write(*,*)
    end if
    assembletime=omp_get_wtime()

    allocate(f(nvars),x(nvars))
    f=0
    x=0

    allocate(A(nvars,nvars))
    A=0
    node_xy_t = transpose(node_xy(:,2:3))
    ele_num_t = transpose(ele_num)
    call assemble_jb(A,f,node_xy_t,ele_num_t,nnode,&
                        nele,nqp,nvars,nu,u_position,v_position,p_position)                
       
    

    write(*,'(a,f11.7,a)')'Initial Assembly time was',omp_get_wtime()-assembletime,' seconds'
    write(22,'(a,f11.7,a)')'Initial Assembly time was',omp_get_wtime()-assembletime,' seconds'
    
    if(.not.quiet)then
        write(*,*)'Finished assembleing coeficient matrix.'
        write(*,*)
        write(*,*)'Determining boundary types...'
        write(*,*)
    end if
    
    call boundary_type(node_xy(:,2:3),onboundary,f,nnode,nvars,u_position,&
                                                         v_position,p_position)
 
    if(.not.quiet)then
        write(*,*)'Fininished determining boundary types.'
        write(*,*)'Calling boundary condition subroutine...'    
    end if


    call boundary_conditions(A,onboundary,nnode,nvars,u_position,&
                                    v_position,p_position,'boundary_type.out')
 
    if(.not.quiet)then
        write(*,*)'Wrote boundary types to: boundary_type.out'
        write(*,*)'Fininished imposing boundary conditions.'
    end if   

    if(sparse)then  
      
        c1=0
        do i=1,nvars
            do j=1,nvars
                if(A(i,j)/=0d0)then
                    c1=c1+1
                end if
            end do
        end do
        nelt=c1

        allocate(ia(nelt),ja(nelt),aa(nelt))        
 
            c1=0
            do i=1,nvars
                do j=1,nvars
                    if(A(i,j)/=0d0)then
                        c1=c1+1
                        ia(c1)=i
                        ja(c1)=j
                        aa(c1)=A(i,j)
                    end if
                end do
            end do

        if(.not.ns)then
            deallocate(A)
        end if
        
        write(*,'(a,i6,a,i10)')'There were',nelt,' nonzero entries out of',nvars**2
        write(*,'(a,f11.7,a)')'That is',100d0-100d0*dble(nvars)/dble(nelt)**2,'% sparse'
        write(22,'(a,i6,a,i10)')'There were',nelt,' nonzero entries out of',nvars**2
        write(22,'(a,f11.7,a)')'That is',100d0-100d0*dble(nvars)/dble(nelt)**2,'% sparse'
         
        write(*,*)
        solvetime=omp_get_wtime()
        !if(.not.quiet)then
            write(*,*)'Laying the SLAP down with DGMRES..' 
            write(*,*)
        !end if
        
        !alpha=1d0
        !call mkl_dcoosv('n',nvars,alpha,'gunfgg',aa,ia,ja,nelt,f,x)
        
        call slapit(aa,ia,ja,f,nelt,nvars,info,quiet)
        x=f
       
 
        write(*,'(a,f11.7,a)')'Slap time was',omp_get_wtime()-solvetime,' seconds'
        write(22,'(a,f11.7,a)')'Slap time was',omp_get_wtime()-solvetime,' seconds'
        write(122,*)nvars,omp_get_wtime()-solvetime

        if(info/=0.and..not.quiet)then
            write(*,*)
            write(*,*)'IERR = ',info
            write(*,*)'      IERR = 0 => All went well.                       '
            write(*,*)'      IERR = 1 => Insufficient storage allocated for   '
            write(*,*)'                  RGWK or IGWK.                        '
            write(*,*)'      IERR = 2 => Routine Dgmres failed reducing the   '
            write(*,*)'                  norm of the current residual on      '
            write(*,*)'                  its last call,                       '
            write(*,*)'                  and so the iteration has stalled.  In'
            write(*,*)'                  this case, X equals the last computed'
            write(*,*)'                  approximation.  The user must either '
            write(*,*)'                  increase MAXL, or choose a different '
            write(*,*)'                  initial guess.                       '
            write(*,*)'      IERR =-1 => Insufficient length for RGWK array.  '
            write(*,*)'                  IGWK(6) contains the required minimum'
            write(*,*)'                  length of the RGWK array.            '
            write(*,*)'      IERR =-2 => Inconsistent ITOL and JPRE values.   '
            write(*,*)'For IERR <= 2, RGWK(1) = RHOL, which is the norm on the'
            write(*,*)'left-hand-side of the relevant stopping test defined   '
            write(*,*)'below associated with the residual for the current     '
            write(*,*)'approximation X(L).                                    '
            write(*,*)
        end if


    else
        lda=nvars
        ldb=nvars
        nrhs=1
        allocate(ipiv(nvars))
        
        solvetime=omp_get_wtime()

        !if(.not.quiet)then
            write(*,*)'Crunching with LAPACK...'
        !end if

        call DGETRF(nvars,nvars,A,lda,ipiv,info)

        if(.not.quiet)then
            write(*,*)'LAPACK DGETRF exit status was',info,'should be 0'
        end if
        
        call DGETRS('N',nvars,nrhs,A,lda,ipiv,f,ldb,info )
        x=f
        if(.not.quiet)then
            write(*,*)'LAPACK DGETRF exit status was',info,'should be 0'
            write(*,*)
        end if
        
        write(*,'(a,f11.7,a)')'Lapack time was',omp_get_wtime()-solvetime,' seconds'
        write(22,'(a,f11.7,a)')'Lapack time was',omp_get_wtime()-solvetime,' seconds'
        write(122,*)nvars,omp_get_wtime()-solvetime

    end if
       


 

    if(ns)then
        !solve the full nonlinear navier stokes equations with the 
        !stokes solution as the initial guess

       
            allocate(r(nvars),delx(nvars))
            if(sparse)then
                deallocate(ia,ja,aa)        
            end if
            ns_it=0

      if(.not.jb)then  
    
            if(sparse)then
                write(*,*)"Nonlinear NS via sparse matrix storage (coordinate format)"
            else
                write(*,*)"Nonlinear NS via coordinate direct dense matrix storage"
            end if
                
        do
                
            write(*,*)'Some of the solution' 
            call print_double_mat(x(1:12),1,12)
            

            call residual_fem ( nnode, node_xy_t, nele, ele_num_t,&
             nqp, u_position, v_position, p_position, nvars, nu, x, r )
            
           ! if(ns_it==0)then
           !     open(102,file='residual')
           !     open(104,file='positions')
           !     open(105,file='ele_nums')
           !     open(106,file='node_nums')
           !     open(107,file='stokes_solution')
           !     write(*,*)'nu=',nu,'nvars=',nvars,'nnode=',nnode,'nele=',nele
           !     do i=1,nvars
           !         write(102,*)r(i)
           !         write(107,*)x(i)    
           !     end do
           !     do i=1,nele
           !         write(105,*)(ele_num_t(j,i),j=1,6)
           !     end do
           !     do i=1,nnode
           !         write(104,*)u_position(i),v_position(i),p_position(i)
           !         write(106,*)(node_xy_t(j,i),j=1,2)
           !     end do
           !     close(102)
           !     close(104)
           !     close(105)
           !     close(106)
           !     close(107)
           ! end if
            

            write(*,*)'Some of the residual'
            call print_double_mat(r(1:12),1,12)
        
            call boundary_adjust_residual(node_xy(:,2:3), onboundary, r, x, nnode, nvars, &
                u_position, v_position, p_position)

            write(*,*)'Residual Norm',sqrt(sum(r**2d0))
            if(sqrt(sum(r**2d0)) < nstol )then
                write(*,*)'Newton-Rhapson Converged!'
                deallocate(A)
                exit
            else if(ns_it >= 5)then
                write(*,*)'REACHED ITERATION LIMIT!!!!!'
                exit
            end if

            ns_it = ns_it + 1  

            A=0
            call jacobian_fem ( nnode, node_xy_t, nele, &
                 ele_num_t, nqp, u_position, v_position, &
                 p_position, nvars, nu, x, A)
            

            call boundary_conditions(A,onboundary,nnode,nvars,u_position,&
                                    v_position,p_position,'boundary_type.out')            
           
            !if(ns_it==1)then
            !    open(101,file='jacobian')   
            !    do i=1,nvars
            !        write(101,'(10000f10.3)')(a(i,j),j=1,nvars)
            !    end do
            !    close(101)
            !end if

            if(sparse)then
                if(ns_it==1)then

                    c1=0
                    do i=1,nvars
                        do j=1,nvars
                            if(A(i,j)/=0d0)then
                                c1=c1+1
                            end if
                        end do
                    end do
                    nelt=c1

                    allocate(ia(nelt),ja(nelt),aa(nelt)) 

                else if(ns_it/=1.and.nelt/=c1)then
                    write(*,*)'FUCK SHIT BITCH'
                end if 

                c1=0
                do i=1,nvars
                    do j=1,nvars
                        if(A(i,j)/=0d0)then
                            c1=c1+1
                            ia(c1)=i
                            ja(c1)=j
                            aa(c1)=A(i,j)
                        end if
                    end do
                end do
 

                write(*,'(a,i5,a,i10,a,i10)')'Iteration',ns_it,' There were',nelt,' nonzero entries out of',nvars**2
                write(*,'(a,f11.7,a)')'That is',100d0-100d0*dble(nvars)/dble(nelt)**2,'% sparse'
                write(22,'(a,i5,a,i10,a,i10)')'Iteration',ns_it,' There were',nelt,' nonzero entries out of',nvars**2
                write(22,'(a,f11.7,a)')'That is',100d0-100d0*dble(nvars)/dble(nelt)**2,'% sparse'
                 
                write(*,*)
                solvetime=omp_get_wtime()
                if(.not.quiet)then
                    write(*,*)'Laying the SLAP down with DGMRES..' 
                end if
                
                !alpha=1d0
                !call mkl_dcoosv('n',nvars,alpha,'gunfgg',aa,ia,ja,nelt,f,x)
                
                delx = -r
                call slapit(aa,ia,ja,delx,nelt,nvars,info,quiet)

                write(*,'(a,i5,a,f11.7,a)')'Iteration',ns_it,' Slap time was',omp_get_wtime()-solvetime,' seconds'
                write(22,'(a,i5,a,f11.7,a)')'Iteration',ns_it,' Slap time was',omp_get_wtime()-solvetime,' seconds'

            else if(.not.sparse)then

                solvetime=omp_get_wtime()

                if(.not.quiet)then
                    write(*,*)'Calling LAPACK to solve system of equations...'
                end if

                call DGETRF(nvars,nvars,A,lda,ipiv,info)

                !if(.not.quiet)then
                    write(*,*)'LAPACK DGETRF exit status was',info,'should be 0'
                !end if

                delx = -r

                call DGETRS('N',nvars,nrhs,A,lda,ipiv,delx,ldb,info )
                
                !if(.not.quiet)then
                    write(*,'(a,f11.7,a)')'LAPACK DGETRF exit status was',info,'should be 0'
                    write(*,*)
                !end if

                write(*,'(a,f11.7,a)')'Lapack time was',omp_get_wtime()-solvetime,'seconds'
                write(22,'(a,f11.7,a)')'Lapack time was',omp_get_wtime()-solvetime,'seconds'           
            end if
 
            !delx = r
            !write(*,*)'Update Norm',sqrt(sum(delx**2d0))

            x = x + delx

            


        end do

        !else if(jb)then
        !        write(*,*)"Nonlinear NS via Banded Storage"
        !        
        !        deallocate(A)
        !        if(sparse)then
        !            allocate(ipiv(nvars))
        !        end if
        !    
        !        call bandwidth ( 6, nele, ele_num_t, &
        !            nnode, p_position, u_position, v_position, ib )

        !        allocate(A(3*ib+1,nvars))
        !        write (*,'(a,i3)' )' Iteration',ns_it
        !      do
        !             call residual_fem_jb (nnode, node_xy_t, nele, &
        !               ele_num_t, nqp, u_position, &
        !               v_position, p_position, nvars, nu, x, r )
        !         
        !         
        !            call boundary_adjust_residual(node_xy(:,2:3), onboundary, r, x, nnode, nvars, &
        !                u_position, v_position, p_position)
 
        !         
        !             write ( *, '(a,g14.7)' ) ' norm FEM residual = ', sqrt(sum(r**2))
        !         
        !             if ( sqrt(sum(r**2)) < nstol) then
        !               write ( *, '(a)' ) 'Newton-Rhapson Converged!.'
        !               exit
        !             end if
        !         
        !             if ( ns_it >= 5 ) then
        !               write ( *, '(a)' ) 'TOO MANY ITERATIONS!!!'
        !               exit
        !             end if
        !         
        !             ns_it = ns_it + 1
        !             write ( *,'(a,i3)' ) ' Iteration',ns_it
        !
        !                !  Compute the finite element jacobian matrix.
        !             call jacobian_fem_jb ( nnode, node_xy_t, nele, &
        !               ele_num_t, nqp, u_position, &
        !               v_position, p_position, nvars, nu, x, &
        !               ib, a )
        !         
        !                !  Adjust the jacobian for boundary conditions.
        !            call jacobian_adjust_dirichlet ( nnode, node_xy_t, onboundary, u_position,&
        !                        v_position,p_position, nvars, ib, a )

        !                !     Factor the jacobian.
        !             call dgb_fa ( nvars, ib, ib, a, ipiv, info )
        !         
        !             if ( info /= 0 ) then
        !               write ( *, '(a)' ) '  The Jacobian matrix is singular.'
        !               stop
        !             end if
        !                
        !                !  Set up and solve the Newton system J * dX = - F(x)
        !             delx = -r
        !         
        !             call dgb_sl ( nvars, ib, ib, a, ipiv, delx, 0 )
        !         
        !             x = x - delx
        !         
        !         
        !             !node_c_del_norm = sqrt ( sum ( node_c_del(1:variable_num)**2 ) )
        !         
        !         
        !      end do
        end if
    end if


    open(14,file='velocity.out')
    do i=1,nnode
        write(14,*)x(u_position(i)),x(v_position(i))
    end do
    if(.not.quiet)then
        write(*,*)'Wrote velocity information to velocity.out'
    end if
    close(14)


    open(14,file='solution.out')
    do i=1,nnode
        write(14,*)x(u_position(i)),x(v_position(i)),x(p_position(i))
    end do
    if(.not.quiet)then
        write(*,*)'Wrote solution information to solution.out'
    end if
    close(14)

    !open(14,file='pressure.out')
    !do i=1,nnode
    !    write(14,*)f(p_position(i))
    !end do
    !if(.not.quiet)then
    !    write(*,*)'Wrote pressure information to pressure.out'
    !    write(*,*)
    !end if
    !close(14)

    write(*,'(a,f11.7,a)')'Finished Simulation, Time was',OMP_GET_WTIME()-time,' seconds.'
    write(22,'(a,f11.7,a)')'Time was',OMP_GET_WTIME()-time,' seconds.'
    close(22)
    close(122)


    stop
    
    
end 

