subroutine getnewnodes(refine,nrefine,file,node_xy,nnode,ele_nums,nele,allnodes,newnodes)
        implicit none
        integer::nnode,nele,nrefine,i
        integer,dimension(nele)::refine
        double precision,dimension(3)::xx,yy
        double precision,dimension(nnode,2)::node_xy
        integer,dimension(nele,3)::ele_nums
        double precision,dimension(nrefine,2)::newnodes
        double precision,dimension(nrefine+nnode,2)::allnodes
        double precision::up,down,left,right,centerh,centerk,radius,xmid,ymid
        integer::point,ninteriornodes,nsidepts,nobspts
        character(len=10)::oshape
        character(len=*)file
        
        !read in information from original mesh
        open(15,file=file)
        read(15,*)up,down,left,right,ninteriornodes,nsidepts,nobspts
        close(15)

        
        point=0

        do i=1,nele
            xx(1)=node_xy(ele_nums(i,1),1)
            xx(2)=node_xy(ele_nums(i,2),1)
            xx(3)=node_xy(ele_nums(i,3),1)
            yy(1)=node_xy(ele_nums(i,1),2)
            yy(2)=node_xy(ele_nums(i,2),2)
            yy(3)=node_xy(ele_nums(i,3),2)
            if(refine(i)==1)then
                point=point+1
                xmid=(1d0/3d0)*(xx(1)+xx(2)+xx(3))
                if(xmid>right.or.xmid<left)then
                    write(*,*)'out of bounds',i
                end if
                ymid=(1d0/3d0)*(yy(1)+yy(2)+yy(3))
                if(ymid>up.or.ymid<down)then
                    write(*,*)'out of bounds',i
                end if
                newnodes(point,1)=xmid
                newnodes(point,2)=ymid
                allnodes(point+nnode,1)=xmid
                allnodes(point+nnode,2)=ymid
                !write(*,*)nnode+point,newnodes(i,:)
            end if 
        end do
        
        !write(*,*)'-------------- is',point,'the same as',nrefine
        !stick old and new nodes together

!$OMP parallel do
        do i=1,nnode+nrefine
            if(i<=nnode)then
                allnodes(i,1)=node_xy(i,1)
                allnodes(i,2)=node_xy(i,2)
            !elseif(i>nnode)then
            !   write(*,*)i,newnodes(i,:)
            !   allnodes(i,1)=newnodes(i,1)
            !   allnodes(i,2)=newnodes(i,2)
            end if  
        end do
!$OMP end parallel do

return
end