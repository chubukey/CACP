      subroutine build_atom_IVC(g0xyz,x_atoms,maxatom_of_one_node,
     &  nnode_interface,id_node_interface,atom_IVC,natoms,cutoff,
     &  maxnode,maxcrd,natoms_node,weight_atom_IVC,scale_coeff,
     &  xc_IVC,atom_at_node)
      implicit none
      real(kind=8)::xc_IVC(maxnode*3),x_ave(3),xn_tmp(3)
      INTEGER:: maxnode,maxcrd
      REAL(8):: g0xyz(maxcrd)
      INTEGER:: nnode_interface                     ! # of node at atom interface
      INTEGER:: id_node_interface(maxnode)         ! 1~nnode_interface: exact node id
      INTEGER:: natoms_node(maxnode)               ! # of atoms belong to node
      INTEGER:: atom_IVC(maxatom_of_one_node,maxnode)  ! 2d array,  (i,j) ith atom id of node j   
      INTEGER:: atom_at_node(maxnode)
      INTEGER:: natoms
      REAL(kind=8):: w_atom_IVC(maxatom_of_one_node,maxnode)
      REAL(kind=8):: x_atoms(3*natoms)
      integer:: maxatom_of_one_node     
      integer:: node_id,atom_id
      INTEGER:: i,j,k
      real(kind=8):: delx,dely,delz,dist,cutsq,cutoff
      real(kind=8):: xa(3),xn(3),xa_r(3)
      integer:: ncount
      real(kind=8):: ave_weight,scale_coeff,xn0(3),atc_coeff
      REAL(kind=8):: weight_atom_IVC(maxatom_of_one_node,maxnode)
      real(kind=8):: mindis
      
      atc_coeff = 1.0/scale_coeff
      cutsq = cutoff*cutoff*atc_coeff*atc_coeff
!      write(*,*) 'nnode_interface is',nnode_interface
      do i = 1,nnode_interface
         mindis = 1e5
         x_ave(1) = 0.0
         x_ave(2) = 0.0
         x_ave(3) = 0.0
         ncount = 0
         node_id = id_node_interface(i)
         if(node_id == 0) then
            write(*,*) 'error:this node bc is undefined'
         endif
!         write(*,*) 'node_id is',node_id
         xn(1) = g0xyz(node_id*3-2)
         xn(2) = g0xyz(node_id*3-1)
         xn(3) = g0xyz(node_id*3)
!         write(*,*)'node ',xn(1),xn(2),xn(3)
         do j = 1,natoms
            xa(1) = x_atoms(j*3-2)
            xa(2) = x_atoms(j*3-1)
            xa(3) = x_atoms(j*3)
            xa_r(1) = xa(1)*atc_coeff
            xa_r(2) = xa(2)*atc_coeff
            xa_r(3) = xa(3)*atc_coeff
            delx = abs(xn(1)-xa_r(1))
            dely = abs(xn(2)-xa_r(2))
            delz = abs(xn(3)-xa_r(3))
            dist = delx*delx+dely*dely+delz*delz
            if(dist<cutsq) then
               ncount = ncount+1
               atom_IVC(ncount,node_id) = j
               if(dist<mindis) then
                  mindis = dist
                  atom_at_node(node_id) = j
               endif
!               write(*,*) 'node_id,xyz;atom_id,xyz',node_id,
!     &         xn(1),xn(2),xn(3),j,xa(1),xa(2),xa(3)
            endif
         enddo
!        write(*,*) 'natoms_node is',ncount,node_id,
!     & atom_at_node(node_id), mindis
         natoms_node(node_id) = ncount
         ave_weight = 1.0/ncount
         do k = 1,ncount
            weight_atom_IVC(k,node_id) = ave_weight
!            write(*,*) 'weight of atom k on node i is',
!     &  atom_IVC(k,node_id),node_id,weight_atom_of
         enddo
         do k = 1,ncount
            atom_id = atom_IVC(k,node_id)
            x_ave(1)=x_ave(1)+x_atoms(atom_id*3-2)*
     &        weight_atom_IVC(k,node_id)
            x_ave(2)=x_ave(2)+x_atoms(atom_id*3-1)*
     &        weight_atom_IVC(k,node_id)
            x_ave(3)=x_ave(3)+x_atoms(atom_id*3)*
     &        weight_atom_IVC(k,node_id)
         enddo
         xc_IVC(node_id*3-2)  = x_ave(1)*atc_coeff
         xc_IVC(node_id*3-1)  = x_ave(2)*atc_coeff 
         xc_IVC(node_id*3  )  = x_ave(3)*atc_coeff 
!         write(*,*)'x_ave,xn',i,node_id,natoms_node(node_id),xn(1),
!     &  xc_IVC(node_id*3-2),xn(2),xc_IVC(node_id*3-1),
!     &  xn(3),xc_IVC(node_id*3)
      enddo
      
      return 
      end
      
      subroutine Atom_Voronoi_Partion(g0xyz,x_atoms,zu,maxatom_node,
     &  nnode_interface,id_node_interface,thick1,thick2,natoms,me,
     &  natoms_IVC,atom_IVC,xc_IVC,weight_atom_IVC,
     &  natoms_SVC,atom_SVC,xc_SVC,weight_atom_SVC,
     &  natoms_FVC,atom_FVC,xc_FVC,weight_atom_FVC,
     &  maxnode,maxcrd,scale_coeff,r1,r2,group_tag)
      implicit none
      integer:: me
      real(kind=8)::xc_IVC(nnode_interface*3)
      real(kind=8)::xc_SVC(nnode_interface*3)
      real(kind=8)::xc_FVC(nnode_interface*3)
      real(kind=8)::x_ave(3),xn_tmp(3)
      INTEGER:: maxnode,maxcrd
      REAL(8):: g0xyz(maxcrd)
      INTEGER:: nnode_interface                     ! # of node at atom interface
      INTEGER:: id_node_interface(nnode_interface)         ! 1~nnode_interface: exact node id
      INTEGER:: atom_IVC(maxatom_node,nnode_interface)  ! 2d array,  (i,j) ith atom id of node j   
      INTEGER:: atom_SVC(maxatom_node,nnode_interface)  ! 2d array,  (i,j) ith atom id of node j   
      INTEGER:: atom_FVC(maxatom_node,nnode_interface)  ! 2d array,  (i,j) ith atom id of node j   
      INTEGER:: atom_at_node(nnode_interface)
      INTEGER:: natoms
      REAL(kind=8):: x_atoms(3*natoms),zu(natoms)
      INTEGER(kind=4):: group_tag(natoms)
      INTEGER(kind=4):: natoms_IVC(nnode_interface)               ! # of atoms belong to node
      INTEGER(kind=4):: natoms_SVC(nnode_interface)               ! # of atoms belong to node
      INTEGER(kind=4):: natoms_FVC(nnode_interface)               ! # of atoms belong to node
      integer:: maxatom_node,id_nearest_node   
      integer:: node_id,atom_id,id_node,n_all_IN
      INTEGER:: i,j,k
      real(kind=8):: delx,dely,delz,dist,dist1,cutsq
      real(kind=8):: thick1,thick2,r1,r2
      real(kind=8):: xa(3),xn(3),xa_r(3)
      real(kind=8):: x0(3)
      integer:: ncount
      integer:: n_all_IVC,n_IVC_check
      integer:: n_all_SVC,n_SVC_check
      integer:: n_all_FVC,n_FVC_check
      real(kind=8):: ave_weight,scale_coeff,xn0(3),atc_coeff
      REAL(kind=8):: weight_atom_IVC(nnode_interface)
      REAL(kind=8):: weight_atom_SVC(nnode_interface)
      REAL(kind=8):: weight_atom_FVC(nnode_interface)
      real(kind=8):: mindis,distcheck
      real(4):: ave_natoms_IVC,ave_natoms_SVC,ave_natoms_FVC
      real(4):: r4array(10)
      atc_coeff = 1.0/scale_coeff
      x0(1) = 0.d0
      x0(2) = 0.d0
!      write(*,*) 'basic check',r1,thick1,natoms,nnode_interface
      n_all_IVC = 0
      n_all_SVC = 0
      n_all_FVC = 0
      n_all_IN = 0

      ave_natoms_IVC = 0.0
      ave_natoms_SVC = 0.0
      ave_natoms_FVC = 0.0

      do i = 1,nnode_interface
         natoms_SVC(i) = 0
      enddo

      do j = 1,natoms
         group_tag(j) = 0
         !first check if they are in interface region,require the center to be (0,0,0)
         xa(1:3) = x_atoms((j*3-2):(j*3))
         xa(3) = zu(j)            ! modified at May_24, should use zu instead of z to construct the atom_of_node relation
         xa_r(1) = xa(1)*atc_coeff
         xa_r(2) = xa(2)*atc_coeff
         xa_r(3) = xa(3)*atc_coeff
         dist = (sqrt(xa(1)*xa(1)+xa(2)*xa(2))-r1)
         dist1 = (sqrt(xa(1)*xa(1)+xa(2)*xa(2))-r2)
        ! ----------- IVC criteria --------------!
         if(abs(dist)<thick1) then
!            write(*,*) 'check 1'
            group_tag(j) = -j  ! -1 means in IVC, +j means in SVC(j),0 means inside
            n_all_IVC = n_all_IVC+1
            id_nearest_node = 0
            mindis = 1.0e5
            do i = 1,nnode_interface
               node_id = id_node_interface(i)
               if(node_id == 0) then
                  write(*,*) 'error:this node bc is undefined in IVC'
               endif
               xn(1:3) = g0xyz((node_id*3-2):(node_id*3))
               delx = xn(1)-xa_r(1)
               dely = xn(2)-xa_r(2)
               delz = xn(3)-xa_r(3)
               dist = delx*delx+dely*dely+delz*delz
!              if(xa(3)>0.1) then
!              write(*,*) 'atom',j,'node',i,node_id,
!     & 'dist is',xa(1:3),xn(1:3),dist,mindis
!              endif
               if(dist<mindis) then
                  mindis = dist
                  id_nearest_node = i
               endif 
            enddo

            id_node = id_nearest_node   ! indexed from 1->n_inter_node
            natoms_IVC(id_node) =  natoms_IVC(id_node)+1
            atom_IVC(natoms_IVC(id_node),id_node) = j
            group_tag(j) = -id_node
!            write(*,*) 'check 1.1'
        ! ----------- SVC criteria --------------!
         elseif(dist>thick1)then
!            write(*,*) 'check 2'
            n_all_SVC = n_all_SVC+1
            id_nearest_node = 0
            mindis = 1.0e5
            do i = 1,nnode_interface
               node_id = id_node_interface(i)
               if(node_id == 0) then
                  write(*,*) 'error:this node bc is undefined in SVC',
     & id_node_interface(i)
               endif
               xn(1:3) = g0xyz((node_id*3-2):(node_id*3))
               delx = xn(1)-xa_r(1)
               dely = xn(2)-xa_r(2)
               delz = xn(3)-xa_r(3)
               dist = delx*delx+dely*dely+delz*delz
               if(dist<mindis) then
                  mindis = dist
                  id_nearest_node = i
               endif 
            enddo
            id_node = id_nearest_node
            natoms_SVC(id_node) =  natoms_SVC(id_node)+1
            atom_SVC(natoms_SVC(id_node),id_node) = j
            group_tag(j) = id_node

!            write(*,*) 'check 2.1'
        ! ----------- FVC criteria --------------!
          elseif(abs(dist1)<thick2)then
!            write(*,*) 'check 3'
            n_all_FVC = n_all_FVC+1
            id_nearest_node = 0
            mindis = 1.0e5
            do i = 1,nnode_interface
               node_id = id_node_interface(i)
               if(node_id == 0) then
                  write(*,*) 'error:this node bc is undefined in FVC'
               endif
               xn(1:3) = g0xyz((node_id*3-2):(node_id*3))
               delx = xn(1)-xa_r(1)
               dely = xn(2)-xa_r(2)
               delz = xn(3)-xa_r(3)
               dist = delx*delx+dely*dely+delz*delz
               if(dist<mindis) then
                  mindis = dist
                  id_nearest_node = i
               endif 
            enddo
            id_node = id_nearest_node
            natoms_FVC(id_node) =  natoms_FVC(id_node)+1
            atom_FVC(natoms_FVC(id_node),id_node) = j
            group_tag(j) = 0
          else 
            n_all_IN = n_all_IN + 1
          endif
          

!            write(*,*) 'check 3.1'
      enddo
!      write(*,*) 'atom_IVC is',n_all_IVC,'atom_SVC is',n_all_SVC,
!     & 'atom_IN',n_all_IN,n_all_IVC+n_all_SVC+n_all_IN
      ! careful here xc_IVC is not calculated yet      
      n_IVC_check = 0
      n_SVC_check = 0
      n_FVC_check = 0
 
      do i = 1,nnode_interface
         node_id = id_node_interface(i)
         xn(1:3) = g0xyz((node_id*3-2):(node_id*3))
         n_IVC_check = n_IVC_check+natoms_IVC(i)
         n_SVC_check = n_SVC_check+natoms_SVC(i)
         n_FVC_check = n_FVC_check+natoms_FVC(i)
         ave_natoms_IVC = ave_natoms_IVC + natoms_IVC(i)/nnode_interface
         ave_natoms_SVC = ave_natoms_SVC + natoms_SVC(i)/nnode_interface
         ave_natoms_FVC = ave_natoms_FVC + natoms_FVC(i)/nnode_interface
         if(me==0) then 
!          write(*,*) 'n IVC,SVC,FVC is',i,node_id,
!     &        natoms_IVC(i),natoms_SVC(i),natoms_FVC(i)
          r4array(1) = i
          r4array(2) = node_id
          r4array(3) = natoms_IVC(i)
          r4array(4) = natoms_SVC(i)
          r4array(5) = natoms_FVC(i)

          write(*,*) 'n IVC,SVC,FVC is',r4array(2:5)
!          write(*,*) 'n IVC,SVC,FVC is',i
!     &        natoms_IVC(i),natoms_SVC(i),natoms_FVC(i)
         endif
!         write(*,*) 'n IVC,SVC is',i,natoms_IVC(i),natoms_SVC(i),xn(1:3)
!         node_id = id_node_interface(i)
         weight_atom_IVC(i) = 1.0/natoms_IVC(i)
         weight_atom_SVC(i) = 1.0/natoms_SVC(i)
         weight_atom_FVC(i) = 1.0/natoms_FVC(i)

      ! ------------ IVC part ------------ !

         x_ave(1) = 0.d0
         x_ave(2) = 0.d0
         x_ave(3) = 0.d0
         do k = 1,natoms_IVC(i)
            atom_id = atom_IVC(k,i);
            x_ave(1)=x_ave(1)+x_atoms(atom_id*3-2)*weight_atom_IVC(i)
            x_ave(2)=x_ave(2)+x_atoms(atom_id*3-1)*weight_atom_IVC(i)
            x_ave(3)=x_ave(3)+x_atoms(atom_id*3-0)*weight_atom_IVC(i)
         enddo          
         xc_IVC(i*3-2) = x_ave(1)*atc_coeff
         xc_IVC(i*3-1) = x_ave(2)*atc_coeff
         xc_IVC(i*3-0) = x_ave(3)*atc_coeff


      ! ------------ SVC part ------------ !
         x_ave(1) = 0.d0
         x_ave(2) = 0.d0
         x_ave(3) = 0.d0
         do k = 1,natoms_SVC(i)
            atom_id = atom_SVC(k,i);
            x_ave(1)=x_ave(1)+x_atoms(atom_id*3-2)*weight_atom_SVC(i)
            x_ave(2)=x_ave(2)+x_atoms(atom_id*3-1)*weight_atom_SVC(i)
            x_ave(3)=x_ave(3)+x_atoms(atom_id*3-0)*weight_atom_SVC(i)
         enddo          
         xc_SVC(i*3-2) = x_ave(1)*atc_coeff
         xc_SVC(i*3-1) = x_ave(2)*atc_coeff
         xc_SVC(i*3-0) = x_ave(3)*atc_coeff


      ! ------------ FVC part ------------ !
         x_ave(1) = 0.d0
         x_ave(2) = 0.d0
         x_ave(3) = 0.d0
         do k = 1,natoms_FVC(i)
            atom_id = atom_FVC(k,i);
            x_ave(1)=x_ave(1)+x_atoms(atom_id*3-2)*weight_atom_FVC(i)
            x_ave(2)=x_ave(2)+x_atoms(atom_id*3-1)*weight_atom_FVC(i)
            x_ave(3)=x_ave(3)+x_atoms(atom_id*3-0)*weight_atom_FVC(i)
         enddo          
         xc_FVC(i*3-2) = x_ave(1)*atc_coeff
         xc_FVC(i*3-1) = x_ave(2)*atc_coeff
         xc_FVC(i*3-0) = x_ave(3)*atc_coeff
      enddo    
      
   
      if(me==0) then
         write(*,*) 'ave_natoms_IVC/SVC/FVC is',
     &   ave_natoms_IVC,ave_natoms_SVC,ave_natoms_FVC
         write(*,*) 'atom_IVC is',n_all_IVC,
     &     'atom_SVC is',n_all_SVC,'atom_FVC is',n_all_FVC
         write(*,*) 'check_IVC is',n_IVC_check,
     &      'check_SVC is',n_SVC_check,'check_FVC is',n_FVC_check
      endif

      return 
      end
      
      
      subroutine get_npadding_atoms(natoms,x_atoms,npadding_atoms,
     & MD_box_hi,MD_box_lo,scale_coeff)
      implicit none
      real(kind=8)::MD_box_hi(3),MD_box_lo(3)
      integer::i,natoms,npadding_atoms
      real(kind=8)::x_atoms(3*natoms)
      real(kind=8)::xa(3),xa_r(3)
      integer:: ncount
      real(kind=8)::scale_coeff,atc_coeff
      
      atc_coeff = 1.0/scale_coeff
      ncount = 0      
    
      do i = 1,natoms
         xa(1) = x_atoms(i*3-2)
         xa(2) = x_atoms(i*3-1)
         xa(3) = x_atoms(i*3)
         xa_r(1) = xa(1)*atc_coeff
         xa_r(2) = xa(2)*atc_coeff
         xa_r(3) = xa(3)*atc_coeff
         if((xa_r(1)>MD_box_hi(1)).OR.(xa_r(1)<MD_box_lo(1)).OR.
     & (xa_r(2)>MD_box_hi(2)).OR.(xa_r(2)<MD_box_lo(2)).OR.
     & (xa_r(3)>MD_box_hi(3)).OR.(xa_r(3)<MD_box_lo(3))) then
            ncount = ncount+1
!            write(*,*) i,xa_r(1),xa_r(2),xa_r(3)
         endif
      enddo

      npadding_atoms = ncount
      return
      end



      subroutine build_padding_atoms_1(natoms,x_atoms,npadding_atoms,
     & MD_box_hi,MD_box_lo,P_r_lo,P_r_hi,
     & padding_atoms,elem_of_atom, maxcrd,g0xyz,ijk,nx,nelx,mnelx,
     & scale_coeff,shape_coeff_atom)

      implicit none
      real(kind=8)::scale_coeff,atc_coeff
      real(kind=8)::MD_box_hi(3),MD_box_lo(3)
      real(kind=8)::P_r_hi(3),P_r_lo(3)
      integer::i,j,natoms,npadding_atoms
      real(kind=8)::x_atoms(3*natoms)
      real(kind=8)::xa(3),xa_r(3)
      integer:: ncount
      integer:: padding_atoms(npadding_atoms)
      integer:: elem_id,nx,nelx,mnelx,maxcrd
      integer:: ijk(mnelx),nodeset(4)
      real(kind=8)::g0xyz(maxcrd),xset(4),yset(4),zset(4)
      integer:: elem_of_atom(npadding_atoms)
      real(kind=8)::shape_coeff_atom(npadding_atoms*3)
      real(kind=8)::shapev(4),e_center(3)
      integer:: num_boundary_elem, boundary_elem(nelx),elem_count
      LOGICAL:: flag_in(3)
      ! target: find padding_atoms(),elem_of_atom(),shape_coeff_atom()
      ncount = 0      
      open(355,file='padding_atoms.dat')
!      write(355,*)'count,atom_id,elem_id,shapev(1),shapev(2),shapev(3)'
! first construct the possible elements contain the padding elements, should be much less than nelements
      elem_count = 0
      do i=1,nelx
         boundary_elem(i) = 0
      enddo

      do i=1,nelx
         e_center(1) = 0.0
         e_center(2) = 0.0
         e_center(3) = 0.0
         nodeset(1)=ijk(i*4-3)
         nodeset(2)=ijk(i*4-2)
         nodeset(3)=ijk(i*4-1)
         nodeset(4)=ijk(i*4)
      
         do j=1,4
            xset(j)=g0xyz(nodeset(j)*3-2)
            yset(j)=g0xyz(nodeset(j)*3-1)
            zset(j)=g0xyz(nodeset(j)*3)
            e_center(1) = e_center(1)+xset(j)/4;
            e_center(2) = e_center(2)+yset(j)/4;
            e_center(3) = e_center(3)+zset(j)/4;
          enddo
          
          flag_in(1)=(e_center(1)>P_r_lo(1)).AND.(e_center(1)<P_r_hi(1))
          flag_in(2)=(e_center(2)>P_r_lo(2)).AND.(e_center(2)<P_r_hi(2))
          flag_in(3)=(e_center(3)>P_r_lo(3)).AND.(e_center(3)<P_r_hi(3))
          if((flag_in(1).AND.flag_in(2).AND.flag_in(3))) then
             elem_count = elem_count + 1
             boundary_elem(elem_count) = i
          endif
      enddo
      num_boundary_elem = elem_count
      write(*,*) 'possible padding elements # is',num_boundary_elem
 
      ! careful, not sure about whether i is atom_id yet, kind of confirm
      do i = 1,natoms
         xa(1) = x_atoms(i*3-2)
         xa(2) = x_atoms(i*3-1)
         xa(3) = x_atoms(i*3)
         xa_r(1) = xa(1)*atc_coeff
         xa_r(2) = xa(2)*atc_coeff
         xa_r(3) = xa(3)*atc_coeff
         if((xa_r(1)>MD_box_hi(1)).OR.(xa_r(1)<MD_box_lo(1)).OR.
     & (xa_r(2)>MD_box_hi(2)).OR.(xa_r(2)<MD_box_lo(2)).OR.
     & (xa_r(3)>MD_box_hi(3)).OR.(xa_r(3)<MD_box_lo(3))) then
            ncount = ncount+1
            padding_atoms(ncount) = i
            call findpointvalue(xa_r,g0xyz,ijk,nx,nelx,elem_id,shapev,
     & i,boundary_elem,num_boundary_elem)
            elem_of_atom(ncount) = elem_id
            shape_coeff_atom(ncount*3-2) = shapev(1)
            shape_coeff_atom(ncount*3-1) = shapev(2)
            shape_coeff_atom(ncount*3)   = shapev(3)
             
            write(355,*) ncount,i,elem_id,shapev(1),shapev(2),shapev(3)
         endif
      enddo
      close(355)

      return
      end

      subroutine build_padding_atoms_2(natoms,x_atoms,npadding_atoms,
     & MD_box_hi,MD_box_lo,P_r_lo,P_r_hi,
     & padding_atoms,elem_of_atom, maxcrd,g0xyz,ijk,nx,nelx,mnelx,
     & scale_coeff,shape_coeff_atom)

      implicit none
      real(kind=8)::scale_coeff,atc_coeff
      real(kind=8)::MD_box_hi(3),MD_box_lo(3)
      real(kind=8)::P_r_hi(3),P_r_lo(3)
      integer::i,j,natoms,npadding_atoms
      real(kind=8)::x_atoms(3*natoms)
      real(kind=8)::xa(3),xa_r(3)
      integer:: ncount
      integer:: padding_atoms(npadding_atoms)
      integer:: elem_id,nx,nelx,mnelx,maxcrd
      integer:: ijk(mnelx),nodeset(4)
      real(kind=8)::g0xyz(maxcrd),xset(4),yset(4),zset(4)
      integer:: elem_of_atom(npadding_atoms)
      real(kind=8)::shape_coeff_atom(npadding_atoms*3)
      real(kind=8)::shapev(4),e_center(3),x_check(3)
      integer:: num_boundary_elem, boundary_elem(nelx),elem_count
      LOGICAL:: flag_in(3)
      ! target: find padding_atoms(),elem_of_atom(),shape_coeff_atom()
      ncount = 0      
      atc_coeff = 1.0/scale_coeff
! first construct the possible elements contain the padding elements, should be much less than nelements
      elem_count = 0
      do i=1,nelx
         boundary_elem(i) = 0
      enddo

      do i=1,nelx
         e_center(1) = 0.0
         e_center(2) = 0.0
         e_center(3) = 0.0
         nodeset(1)=ijk(i*4-3)
         nodeset(2)=ijk(i*4-2)
         nodeset(3)=ijk(i*4-1)
         nodeset(4)=ijk(i*4)
      
         do j=1,4
            xset(j)=g0xyz(nodeset(j)*3-2)
            yset(j)=g0xyz(nodeset(j)*3-1)
            zset(j)=g0xyz(nodeset(j)*3)
            e_center(1) = e_center(1)+xset(j)/4;
            e_center(2) = e_center(2)+yset(j)/4;
            e_center(3) = e_center(3)+zset(j)/4;
          enddo
          
          flag_in(1)=(e_center(1)>P_r_lo(1)).AND.(e_center(1)<P_r_hi(1))
          flag_in(2)=(e_center(2)>P_r_lo(2)).AND.(e_center(2)<P_r_hi(2))
          flag_in(3)=(e_center(3)>P_r_lo(3)).AND.(e_center(3)<P_r_hi(3))
          if((flag_in(1).AND.flag_in(2).AND.flag_in(3))) then
             elem_count = elem_count + 1
             boundary_elem(elem_count) = i
          endif
      enddo
      num_boundary_elem = elem_count
      write(*,*) 'possible padding elements # is',num_boundary_elem
 
      ! careful, not sure about whether i is atom_id yet, kind of confirm
      do i = 1,natoms
         x_check(1) = 0.0d0
         x_check(2) = 0.0d0
         x_check(3) = 0.0d0
         xa(1) = x_atoms(i*3-2)
         xa(2) = x_atoms(i*3-1)
         xa(3) = x_atoms(i*3)
         xa_r(1) = xa(1)*scale_coeff
         xa_r(2) = xa(2)*scale_coeff
         xa_r(3) = xa(3)*scale_coeff
         if((xa_r(1)>MD_box_hi(1)).OR.(xa_r(1)<MD_box_lo(1)).OR.
     & (xa_r(2)>MD_box_hi(2)).OR.(xa_r(2)<MD_box_lo(2)).OR.
     & (xa_r(3)>MD_box_hi(3)).OR.(xa_r(3)<MD_box_lo(3))) then
            ncount = ncount+1
            padding_atoms(ncount) = i
            call findpointvalue(xa_r,g0xyz,ijk,nx,nelx,elem_id,shapev,
     & i,boundary_elem,num_boundary_elem)
            elem_of_atom(ncount) = elem_id
            shape_coeff_atom(ncount*3-2) = shapev(1)
            shape_coeff_atom(ncount*3-1) = shapev(2)
            shape_coeff_atom(ncount*3)   = shapev(3)
            
         endif
      enddo

      return
      end

      subroutine read_padding_atoms(npadding_atoms,padding_atoms,
     &            elem_of_atom,shape_coeff_atom)
      implicit none 
      integer:: npadding_atoms
      integer:: padding_atoms(npadding_atoms)
      integer:: elem_of_atom(npadding_atoms)
      real(kind=8)::shape_coeff_atom(npadding_atoms*3)
      integer:: file_status,i,atom_id,ncount
      
      open(unit=357,file='new_padding_atoms.dat')
      read(357,*) file_status
      write(*,*) 'npadding_atoms is ',npadding_atoms
      do i = 1,npadding_atoms
         read(357,*) ncount,atom_id,elem_of_atom(i),
     & shape_coeff_atom(i*3-2),shape_coeff_atom(i*3-1),
     & shape_coeff_atom(i*3)
         if((i-ncount).ne.0) then
            write(*,*) 'wrong reading'
         endif
         padding_atoms(ncount) = atom_id
      enddo
      write(*,*) 'id of padding atom 20k is ',padding_atoms(20000)
      close(357)

      return
      end


      subroutine ATC_update(u,index_inter,natoms,
     &  avedisp_x,avedisp_y,avedisp_z,
     &  nx,nodeid2interid,
     &  pbc_box,pbc_strain,iflag_topbot,
     &  maxatom_of_one_node,nnode_inter,internode_pbc,scale_coeff,
     &  natoms_IVC,atom_IVC,w_atom_IVC,xc_IVC,
     &  natoms_FVC,atom_FVC,w_atom_FVC,xc_FVC)
      implicit none
      integer::iproc
      real(kind=8):: u(3),xn(3),coeff,scale_coeff,tmp_coeff,xa_tmp(3)
      integer:: natoms,index_inter
      integer:: nx,nodeid2interid(nx),g0xyz(nx*3)
      integer(kind=4):: natoms_IVC(nnode_inter)
      integer(kind=4):: natoms_FVC(nnode_inter)
      INTEGER:: atom_IVC(maxatom_of_one_node,nnode_inter)  ! 2d array,  (i,j) ith atom id of node j   
      INTEGER:: atom_FVC(maxatom_of_one_node,nnode_inter)  ! 2d array,  (i,j) ith atom id of node j   
      REAL(kind=8):: w_atom_IVC(nnode_inter)
      REAL(kind=8):: w_atom_FVC(nnode_inter)
      REAL(kind=8)::xc_IVC(nnode_inter*3)
      REAL(kind=8)::xc_FVC(nnode_inter*3)
      REAL(kind=8)::avedisp_x(natoms)
      REAL(kind=8)::avedisp_y(natoms)
      REAL(kind=8)::avedisp_z(natoms)
      integer:: maxatom_of_one_node,nnode_inter
      integer::internode_pbc(nnode_inter),id_pbc,index_pbc
      integer:: i,j,atom_id,test11
      integer:: iflag_FVCcalc,iflag_topbot
      real(kind=8)::xa_ave(3),atc_coeff,xa_ave_pbc(3)
      real(kind=8)::pbc_box(3),pbc_strain(3),val_topbot
      common/processorinfo/iproc
!      open(283,position = 'Append',file='ATC_update.dat')
!      if(index_inter == 2)then
!        write(*,*) 'pbc box is',pbc_box(1:3)
!        write(*,*) 'pbc strain is',pbc_strain(1:3)
!      endif
!      write(*,*) 'ATC check natoms',index_inter,natoms_IVC(index_inter),
!     & w_atom_IVC(index_inter),scale_coeff
      atc_coeff = 1.0/scale_coeff
!      write(*,*) 'atc_coeff is',atc_coeff
      xa_ave(1:3) = 0.0
      xa_ave_pbc(1:3) = 0.0
!      xa_tmp(1) = xc_IVC(index_inter*3-2)
!      xa_tmp(2) = xc_IVC(index_inter*3-1)
!      xa_tmp(3) = xc_IVC(index_inter*3)
!      write(*,*) 'xn_tmp(1) is',xa_tmp(1),xa_tmp(2),xa_tmp(3)
!      if((index_inter == 2).or.(index_inter == 193))then
!         write(*,*) 'index:',index_inter,'flag',iflag_topbot
!      endif
      val_topbot = 0.5*pbc_strain(3)*pbc_box(3)
      do j=1,natoms_IVC(index_inter)
         atom_id = atom_IVC(j,index_inter)
         if(atom_id == 0) then
             write(*,*)'atom not found IVC',j,index_inter
         endif
         tmp_coeff = w_atom_IVC(index_inter)*atc_coeff
         xa_ave(1)=xa_ave(1)+avedisp_x(atom_id)*tmp_coeff
         xa_ave(2)=xa_ave(2)+avedisp_y(atom_id)*tmp_coeff
         xa_ave(3)=xa_ave(3)+avedisp_z(atom_id)*tmp_coeff
      enddo
      xc_IVC(index_inter*3-2) = xa_ave(1)
      xc_IVC(index_inter*3-1) = xa_ave(2)
      xc_IVC(index_inter*3-0) = xa_ave(3)
      coeff = 1

      id_pbc = internode_pbc(index_inter)
      if(id_pbc.ne.0)then
        index_pbc =  nodeid2interid(id_pbc)
        do j=1,natoms_IVC(index_pbc)
          atom_id = atom_IVC(j,index_pbc)
          if(atom_id == 0) then
              write(*,*)'atom not found IVC',j,index_pbc
          endif
          tmp_coeff = w_atom_IVC(index_pbc)*atc_coeff
          xa_ave_pbc(1)=xa_ave_pbc(1)+avedisp_x(atom_id)*tmp_coeff
          xa_ave_pbc(2)=xa_ave_pbc(2)+avedisp_y(atom_id)*tmp_coeff
          xa_ave_pbc(3)=xa_ave_pbc(3)+avedisp_z(atom_id)*tmp_coeff
        enddo
        do i = 1,3
         u(i) = (xa_ave(i)*natoms_IVC(index_inter)
     &          +xa_ave_pbc(i)*natoms_IVC(index_pbc))
     &          /(natoms_IVC(index_pbc)+natoms_IVC(index_inter))
        enddo
        if(iflag_topbot == 1)then ! top surf       
        ! it is wrong to assign same disp z value for all top surf
        ! nodes, there is difference of u_z, the only same part is the
        ! difference from top surf and bot surf, what is the correct way
        ! to implement it

          u(3) = u(3)+val_topbot
        elseif(iflag_topbot == -1)then ! bot surf
          u(3) = u(3)-val_topbot
        endif
      else
         u(1) = (xa_ave(1))
         u(2) = (xa_ave(2))
         u(3) = (xa_ave(3))
      endif
      
      iflag_FVCcalc = 0
      if(iflag_FVCcalc == 1)then
        do j=1,natoms_FVC(index_inter)
           atom_id = atom_FVC(j,index_inter)
           if(atom_id == 0) then
               write(*,*)'atom not found FVC',j,index_inter
           endif
         ! here it is not u, but should be real coord
           tmp_coeff = w_atom_FVC(index_inter)*atc_coeff
           xa_ave(1)=xa_ave(1)+avedisp_x(atom_id)*tmp_coeff
           xa_ave(2)=xa_ave(2)+avedisp_y(atom_id)*tmp_coeff
           xa_ave(3)=xa_ave(3)+avedisp_z(atom_id)*tmp_coeff
        enddo
      
        xc_FVC(index_inter*3-2) = xa_ave(1)
        xc_FVC(index_inter*3-1) = xa_ave(2)
        xc_FVC(index_inter*3-0) = xa_ave(3)
      endif
!      write(*,*) 'xa compare',xa_ave(1),xa_tmp(1),xa_ave(2),
!     & xa_tmp(2),xa_ave(3),xa_tmp(3)
!      write(*,*) 'u123',xc_IVC(node_id*3-2)-xn_tmp(1)
!     write(*,*) 'u value by ATC',index_inter,u(1),u(2),u(3)
!      write(283,*) index_inter,u(1),u(2),u(3)
!      close(283)
      return
      end
     

      subroutine CTA_u_update(natoms,x_atoms,elem_of_atom,
     & npadding_atoms,padding_atoms,gxyz,shape_coeff_atom,
     & scale_coeff,ijk,mnelx,maxcrd)
      implicit none
     
      real(kind=8)::scale_coeff,atc_coeff
      real(kind=8)::Inter_node_hi(3),Inter_node_lo(3)
      integer::i,j,natoms,npadding_atoms
      real(kind=8)::x_atoms(3*natoms)
      real(kind=8)::xa(3),xa_ori(3),xa_new(3)
      integer:: ncount
      integer:: padding_atoms(npadding_atoms)
      integer:: atom_id,elem_id,nx,nelx,mnelx,maxcrd
      integer:: ijk(mnelx),node_set(4)
      real(kind=8)::gxyz(maxcrd)
      integer:: elem_of_atom(npadding_atoms)
      real(kind=8)::shape_coeff_atom(npadding_atoms*3)
      real(kind=8)::shapev(4)
      atc_coeff = 1.0/scale_coeff
      do i = 1,npadding_atoms
         xa_new(1) = 0.0
         xa_new(2) = 0.0
         xa_new(3) = 0.0
         atom_id =  padding_atoms(i)
         elem_id = elem_of_atom(i)

         xa_ori(1) = x_atoms(atom_id*3-2)*atc_coeff
         xa_ori(2) = x_atoms(atom_id*3-1)*atc_coeff
         xa_ori(3) = x_atoms(atom_id*3)*atc_coeff
         
         node_set(1) = ijk(elem_id*4-3)
         node_set(2) = ijk(elem_id*4-2)
         node_set(3) = ijk(elem_id*4-1)
         node_set(4) = ijk(elem_id*4)

         shapev(1)  = shape_coeff_atom(i*3-2)
         shapev(2)  = shape_coeff_atom(i*3-1)
         shapev(3)  = shape_coeff_atom(i*3)
         shapev(4)  = 1-shapev(1)-shapev(2)-shapev(3)
         if(shapev(4)< -1e-10) then
            write(*,*) 'ERROR: wrong shape_coeff value'
         endif


         do j = 1,4
            xa_new(1) = xa_new(1)+shapev(j)*gxyz(node_set(j)*3-2)
            xa_new(2) = xa_new(2)+shapev(j)*gxyz(node_set(j)*3-1)
            xa_new(3) = xa_new(3)+shapev(j)*gxyz(node_set(j)*3)
         enddo
         if(i == 4000) then
            write(*,*)'compare xa new and old',xa_new(1)-xa_ori(1),
     &      xa_new(1),xa_ori(1)
         endif
!        write(*,*) 'shapev is',shapev(1)+shapev(2)+shapev(3)+shapev(4),
!     &  shapev(1),shapev(2),shapev(3),shapev(4),atom_id
         
         x_atoms(atom_id*3-2) = xa_new(1)/atc_coeff
         x_atoms(atom_id*3-1) = xa_new(2)/atc_coeff
         x_atoms(atom_id*3)   = xa_new(3)/atc_coeff
         
      enddo
      return
      end


      subroutine CTA_f_update(f_react,fa_ext,natoms,mdofx,
     &    maxatom_of_one_node,maxnode,w_atom_IVC,natoms_IVC,
     &    maxcrd,nnode_interface,id_node_interface,atom_IVC,
     &    x_atoms,r,thickness,natoms_SVC,atom_SVC)
      implicit none
     
      real(8):: coeff,scale_coeff,tmp_coeff
      integer:: natoms,node_id
      INTEGER:: nnode_interface                     
      INTEGER:: id_node_interface(nnode_interface)         ! 1~nnode_interface: exact node id
      INTEGER(kind=4):: natoms_IVC(nnode_interface)               ! # of atoms belong to node
      INTEGER(kind=4):: natoms_SVC(nnode_interface)               ! # of atoms belong to node
      INTEGER:: atom_IVC(maxatom_of_one_node,nnode_interface)  ! 2d array,  (i,j) ith atom id of node j   
      INTEGER:: atom_SVC(maxatom_of_one_node,nnode_interface)  ! 2d array,  (i,j) ith atom id of node j   
      REAL(8):: w_atom_IVC(nnode_interface)
      integer:: maxatom_of_one_node,maxnode      
      integer:: i,j
      real(8):: xa_ave(3)
      real(8):: f_react(mdofx),f_node(3)
      integer:: atom_id,elem_id,mdofx,maxcrd
      real(8):: fa_ext(natoms*3),maxfa_ext
      real(8):: unit_coeff,weight
      real(8):: x_atoms(3*natoms),r,thickness,xtmp(3),dist

      ! continuum nodal force unit: MPa * nm^2 -> 1.0e-12 J/m
      ! MD       atomic force unit: eV  / Ang  -> 1.6021766 e-9 J/m
      unit_coeff = 1.0/1602.1766

!      write(*,*) 'CTA_f_update check ',nnode_interface,mdofx
      maxfa_ext = 0.d0
      do i = 1,natoms*3
         maxfa_ext = max(maxfa_ext,abs(fa_ext(i)))
      enddo
!      write(*,*) ' in CTA_f maxfa_ext is',maxfa_ext
      do i = 1,nnode_interface
         node_id = id_node_interface(i)
         f_node(1:3) = f_react((node_id*3-2):(node_id*3))
!         weight = 1.0/(natoms_IVC(i)+natoms_SVC(i))
         weight = 1.0/natoms_IVC(i)
!         write(*,*) 'natoms_IVC SVC is',i,natoms_IVC(i),natoms_SVC(i),
!     & natoms_IVC(i)+natoms_SVC(i),weight
         do j = 1,natoms_IVC(i)
            atom_id = atom_IVC(j,i)
!            fa_ext(atom_id*3-2) = f_node(1)*w_atom_IVC(i)*unit_coeff
!            fa_ext(atom_id*3-1) = f_node(2)*w_atom_IVC(i)*unit_coeff
!            fa_ext(atom_id*3)   = f_node(3)*w_atom_IVC(i)*unit_coeff
            fa_ext(atom_id*3-2) = f_node(1)*weight*unit_coeff
            fa_ext(atom_id*3-1) = f_node(2)*weight*unit_coeff
            fa_ext(atom_id*3-0) = f_node(3)*weight*unit_coeff*0
!         write(*,*) fa_ext(atom_id*3-2),fa_ext(atom_id*3-1),
!     & fa_ext(atom_id*3),f_node(1:3),weight
         enddo
!         do j = 1,natoms_SVC(i)
!            atom_id = atom_SVC(j,i)
!            fa_ext(atom_id*3-2) = f_node(1)*weight*unit_coeff
!            fa_ext(atom_id*3-1) = f_node(2)*weight*unit_coeff
!            fa_ext(atom_id*3-0) = f_node(3)*weight*unit_coeff
!         write(*,*) fa_ext(atom_id*3-2),fa_ext(atom_id*3-1),
!     & fa_ext(atom_id*3),f_node(1:3),weight
!         enddo
            
      enddo
      
      return
      end


      subroutine findpointvalue(ppos,gxyz,ijk,
     & nx,nelx,num_elem,shapev,atom_id,
     & boundary_elem,num_boundary_elem)
      implicit none
      integer::atom_id
      real(8):: gxyz(3*nx),ppos(3),x_check(3)
      integer:: ijk(4*nelx),ijknel(4) 
      real(8):: xset(4),yset(4),zset(4)
      integer::nx,nelx,ii,i,j
      integer::nodeset(4),n1,n2,n3,n4
      real(8):: shapfv(3),shapev(4)
      real(8)::tmp
      integer::nelem_count,num_elem,num_boundary_elem
      integer::boundary_elem(num_boundary_elem)

!   first find which element this point belong to
!   1, find the shape function value of this point in every element
!   2, the element that has all 4 shape function value between 0 and 1 is the element this node belong to 
 
      nelem_count=0          ! incase more than one eleement this point belongs to real8
      num_elem = 0
!      write(*,*) 'num_boundary_elem is',num_boundary_elem
!      do ii=1,num_boundary_elem
!      i = boundary_elem(i)
      do i = 1,nelx     
      nodeset(1)=ijk(i*4-3)
      nodeset(2)=ijk(i*4-2)
      nodeset(3)=ijk(i*4-1)
      nodeset(4)=ijk(i*4)
      
      do j=1,4
      xset(j)=gxyz(nodeset(j)*3-2)
      yset(j)=gxyz(nodeset(j)*3-1)
      zset(j)=gxyz(nodeset(j)*3)
      enddo
!     write(*,*)'xset is',j, xset
!     write(*,*)'yset is',j, yset
!     write(*,*)'zset is',j, zset
      
       
      call findshapefunctionvalue(xset,yset,zset,ppos,shapfv)
!   2
      
      tmp=1.0-shapfv(1)-shapfv(2)-shapfv(3)
      shapev(1)=shapfv(1)
      shapev(2)=shapfv(2)
      shapev(3)=shapfv(3)
      shapev(4)=tmp
      if (((shapfv(2).le.1).and.(shapfv(2).ge.0)).AND.
     & ((shapfv(3).le.1).and.(shapfv(3).ge.0)).AND.
     & ((shapfv(1).le.1).and.(shapfv(1).ge.0)).AND.
     & ((tmp.le.1).and.(tmp.ge.0))) then
      x_check(1) = shapfv(1)*xset(1)+shapfv(2)*xset(2)+shapfv(3)*
     & xset(3)+tmp*xset(4)
      x_check(2) = shapfv(1)*yset(1)+shapfv(2)*yset(2)+shapfv(3)*
     & yset(3)+tmp*yset(4)
      x_check(3) = shapfv(1)*zset(1)+shapfv(2)*zset(2)+shapfv(3)*
     & zset(3)+tmp*zset(4)
       if(x_check(1)-ppos(1)>1e-8)then
          write(*,*)'atom wrong ',x_check(1),ppos(1)
       endif
!      write(*,*)'atom pos',ppos(1),ppos(2),ppos(3),
!     &   shapfv(1),shapfv(2),shapfv(3)
         num_elem = i     
         nelem_count=nelem_count+1
      
!   four nodes of this element      
!   the value at this point
         exit
      endif
      
      
      enddo
      
      if (nelem_count<1) then
         write(*,*)'error, cant find the element for this point',
     & ppos(1),ppos(2),ppos(3),shapev(1),shapev(2),shapev(3),shapev(4)
      else if(nelem_count == 1) then
!         write(*,*)'atom',atom_id,'is in elem:',num_elem
      else
         write(*,*)'error, atom in more than 1 elem'
      endif
      
      return
      end 

      
      
      subroutine findshapefunctionvalue(xset,yset,zset,ppos,shapfv)
      
      implicit none
      real(8):: ppos(3),xset(4), yset(4), zset(4)
      real(8):: shapfv(3),det_A
      real(8)::A(3,3),B(3),A_inv(3,3)
      integer:: i,j
      
      B(1)=ppos(1)-xset(4)
      B(2)=ppos(2)-yset(4)
      B(3)=ppos(3)-zset(4)
            A(1,1)=xset(1)-xset(4)
            A(1,2)=xset(2)-xset(4)
            A(1,3)=xset(3)-xset(4)
            A(2,1)=yset(1)-yset(4)
            A(2,2)=yset(2)-yset(4)
            A(2,3)=yset(3)-yset(4)
            A(3,1)=zset(1)-zset(4)
            A(3,2)=zset(2)-zset(4)
            A(3,3)=zset(3)-zset(4)
     
      
      do i=1,3
      do j=1,3
      A_inv(i,j)=A(i,j) 
      enddo
      enddo
      
      call matinv3(A,A_inv,det_A)
      if(det_A.lt.1.0d-16)then
      write(*,*)'det_A is very small',det_A,'\n'
      endif
      
      call matmult(A_inv,B,shapfv,3,3,1)
      return
      end




       
       
       subroutine matmult(aa,bb,cc,n,m,l)
      
c      include 'aba_param.inc'
       implicit double precision(a-h,o-z)

c   this subroutine returns martix [cc] = matrix [aa] * matrix [bb]
c   n=no. of rows of [aa]      =no. of rows of [cc]
c   m=no. of columns of [aa]   =no. of rows of [bb]
c   l=no. of columns of [bb]   =no. of columns of [cc]

       dimension aa(n,m),bb(m,l),cc(n,l)

       do 10 i=1,n
       do 10 j=1,l
         cc(i,j)=0.0
       do 10 k=1,m
   10    cc(i,j)=cc(i,j)+aa(i,k)*bb(k,j)
       return 
       end
    
      subroutine f_dist_calc(x_atoms,nnode_interface,id_node_interface,
     &  neighlist,n_neigh,natoms,maxnode,atom_at_node,dist)
      implicit none
      real(kind=8)::x_atoms(natoms*3)
      integer:: natoms,maxnode,nnode_interface
      integer:: id_node_interface(maxnode),atom_at_node(maxnode)
      integer:: n_neigh(natoms),neighlist(12,natoms)
       
      integer::i,j,k,node_id,atom_id,neigh_id
      real(kind=8)::x_center(3),x_neigh(3),dist_tmp(3)
      real(kind=8)::dist(12,3,nnode_interface)

!      write(*,*) 'nnode_interface is',nnode_interface
      do i = 1,nnode_interface
         node_id = id_node_interface(i)
         atom_id = atom_at_node(node_id)
         x_center(1) = x_atoms(atom_id*3-2)
         x_center(2) = x_atoms(atom_id*3-1)
         x_center(3) = x_atoms(atom_id*3-0)
    
         do j = 1,12
            dist(j,1,i) = 0.0d0
            dist(j,2,i) = 0.0d0
            dist(j,3,i) = 0.0d0
         enddo
!        write(*,*) 'dist check', node_id,atom_id,n_neigh(atom_id)
         do j = 1,n_neigh(atom_id)
            neigh_id = neighlist(j,atom_id)
            x_neigh(1) = x_atoms(neigh_id*3-2)
            x_neigh(2) = x_atoms(neigh_id*3-1)
            x_neigh(3) = x_atoms(neigh_id*3-0)
             
            dist_tmp(1) = x_neigh(1) - x_center(1)
            dist_tmp(2) = x_neigh(2) - x_center(2)
            dist_tmp(3) = x_neigh(3) - x_center(3)
            dist(j,1,i) = dist_tmp(1)
            dist(j,2,i) = dist_tmp(2)
            dist(j,3,i) = dist_tmp(3)
!            write(*,*) 'dist value',dist_tmp(1),dist_tmp(2),dist_tmp(3)
         enddo

      enddo
!        write(*,*) 'dist value',dist(1,1,1),dist(1,2,1),dist(1,3,1)
      return
      end
 

      subroutine  node_elem_connec_new(ijk,n_elem_node,elem_of_node,
     & maxnode,maxnodelem,iproc,nx,nelx,node,mnelx)
      implicit none
      integer:: elem_of_node(maxnodelem,maxnode),n_elem_node(maxnode)
      integer:: ijk(mnelx)
      integer:: maxnode,maxnodelem,nx,nelx,node,mnelx
      integer:: ncount,i,j,k
      integer:: iproc
!      write(*,*)'nx,nelx,node,mnelx',nx,nelx,node,mnelx,iproc
      do i=1,nx
         ncount = 0
         do j=1,nelx
            do k=1,node
               if(ijk((j-1)*node+k).eq.i)then
                  ncount = ncount+1
                  elem_of_node(ncount,i) = j
!                 write(*,*) 'conn', ncount,ijk((j-1)*node+k),i
!                  write(*,*) 'conn', ncount,ijk((j-1)*node+k),i
               endif
            enddo
         enddo
         n_elem_node(i) = ncount
!         if(iproc == 0) then
!         write(*,*) 'ncount', ncount,n_elem_node(i),i
!         endif
      enddo 
      end

      subroutine  ATC_modify_from_DG_F(DG_FA,DG_FC,nnode_interface,
     & id_node_interface,maxnode,u_CTA,refvector,iproc)
      implicit none
      real(kind=8)::del_DG_F(3,3),DG_FA_tmp(3,3),DG_FC_tmp(3,3)
      integer:: i,j,k,node_id,iproc
      integer:: nnode_interface, maxnode
      integer:: id_node_interface(maxnode)
      real(kind=8)::coords(3)
      real(kind=8)::u_CTA(3,maxnode)
      real(kind=8)::DG_FA(3,3,nnode_interface)
      real(kind=8)::DG_FC(3,3,nnode_interface)
      real(kind=8)::refvector(3,nnode_interface) 
      real(kind=8),parameter::kk=2.0d0
      do i = 1,nnode_interface
         node_id = id_node_interface(i)
         DG_FA_tmp(:,:) = DG_FA(:,:,i)
         DG_FC_tmp(:,:) = DG_FC(:,:,i)
         del_DG_F(:,:) = DG_FC(:,:,i) - DG_FA(:,:,i)
         do j = 1,3
            do k = 1,3
       u_CTA(j,node_id)=kk*del_DG_F(j,k)*refvector(k,i)+u_CTA(j,node_id)
            enddo
         enddo
      enddo

      return
      end
 
      subroutine padding_elem_check(natoms,npadding_atoms,padding_atoms,
     & elem_of_atom,ijk,maxnode,mnelx,nnode_interface,id_node_interface)
      implicit none
      integer:: natoms,npadding_atoms,maxnode,mnelx,nnode_interface
      integer:: elem_of_atom(npadding_atoms)
      integer:: i,j,k,node_id,elem_id,atom_id
      integer:: id_node_interface(maxnode)
      integer:: ijk(mnelx),nodeset(4)
      integer:: padding_atoms(npadding_atoms)
      do i = 1,npadding_atoms
         atom_id = padding_atoms(i)
         elem_id = elem_of_atom(atom_id)
         
         nodeset(1)=ijk(i*4-3)
         nodeset(2)=ijk(i*4-2)
         nodeset(3)=ijk(i*4-1)
         nodeset(4)=ijk(i*4)

         do j = 1,nnode_interface
            node_id = id_node_interface(j)
            if((node_id==nodeset(1)).OR.(node_id==nodeset(2)).OR.
     & (node_id==nodeset(3)).OR.(node_id==nodeset(4))) then
               write(*,*)'warning: ',elem_id,'contains interface_node ',
     & node_id
            exit
            endif
         enddo
      enddo
      return
      end

      subroutine find_inter_node_index(nnode_interface,
     & id_node_interface,nodeid2interid,nx)
      implicit none
      integer::nodeid2interid(nx)
      integer::nnode_interface
      integer::id_node_interface(nnode_interface)
      integer::nx,i,j,id_node
      
      do j = 1,nx
         nodeid2interid(j) = 0
      enddo

      do i = 1,nnode_interface
         id_node = id_node_interface(i)
         nodeid2interid(id_node) = i
      enddo

      return
      end

      subroutine  Surf_f_check(natoms,avedisp_x,avedisp_y,   !node1 is the current node id
     &  avedisp_z,atom_IVC,maxatom_of_one_node,maxnode,nnode_interface,
     &  w_atom_IVC,natoms_IVC,scale_coeff,xc_IVC,iproc)
      implicit none
      integer::iproc
      real(kind=8):: u(3),xn(3),coeff,scale_coeff,tmp_coeff,xa_tmp(3)
      integer:: natoms,index_inter
      integer:: nnode_interface
      INTEGER:: atom_IVC(maxatom_of_one_node,maxnode)  ! 2d array,  (i,j) ith atom id of node j   
      REAL(kind=8):: w_atom_IVC(maxnode)
      REAL(kind=8)::xc_IVC(maxnode*3)
      REAL(kind=8)::avedisp_x(natoms)
      REAL(kind=8)::avedisp_y(natoms)
      REAL(kind=8)::avedisp_z(natoms)
      REAL(kind=8)::aveall_x,aveall_y,aveall_z
      integer:: maxatom_of_one_node,maxnode      
      integer:: i,j,atom_id,test11
      real(kind=8)::xa_ave(3),atc_coeff
      integer(kind=4)::natoms_IVC(maxnode)
      
      open(unit = 137,file = 'surface_check')
      write(137,*) 'u surf check, should be close to 0 for all values'
      write(137,*) 'node_inter #,	u(1),	u(2),	u(3)'
      atc_coeff = 1.0/scale_coeff
      aveall_x = 0.0
      aveall_y = 0.0
      aveall_z = 0.0
      do index_inter = 1,nnode_interface
         xa_ave(1) = 0.0
         xa_ave(2) = 0.0
         xa_ave(3) = 0.0 
         do j=1,natoms_IVC(index_inter)
         
            atom_id = atom_IVC(j,index_inter)
            if(atom_id == 0) then
               write(*,*)'atom not found',j,index_inter
            endif
            ! here it is not u, but should be real coord
          
            tmp_coeff = w_atom_IVC(index_inter)*atc_coeff
            xa_ave(1)=xa_ave(1)+avedisp_x(atom_id)*tmp_coeff
            xa_ave(2)=xa_ave(2)+avedisp_y(atom_id)*tmp_coeff
            xa_ave(3)=xa_ave(3)+avedisp_z(atom_id)*tmp_coeff
         enddo

         xc_IVC(index_inter*3-2) = xa_ave(1)
         xc_IVC(index_inter*3-1) = xa_ave(2)
         xc_IVC(index_inter*3-0) = xa_ave(3)
         coeff = 1
         u(1) = (xa_ave(1))
         u(2) = (xa_ave(2))
         u(3) = (xa_ave(3))
         write(137,*) 'u surf check by ATC',index_inter,u(1),u(2),u(3)
         aveall_x = aveall_x + u(1)*u(1)
         aveall_y = aveall_y + u(2)*u(2)
         aveall_z = aveall_z + u(3)*u(3)
!         write(283,*) index_inter,u(1),u(2),u(3)
!         close(283)
      enddo
      aveall_x = sqrt(aveall_x)
      aveall_y = sqrt(aveall_y)
      aveall_z = sqrt(aveall_z)
      write(137,*) 'sqrt(sum(avedisp^2)) x  y  z'
      write(137,*) aveall_x, aveall_y, aveall_z
      close(137)
      return
      end
    
      subroutine update_MD_boundary(FEM_zlo,FEM_zhi,g0xyz,gpxyz,nx,me)
      implicit none
      integer::i,nzhi,nzlo,nx
      integer::me
      real(kind=8)::zlo,zhi,coeff,FEM_zlo,FEM_zhi
      real(kind=8)::g0xyz(nx*3), gpxyz(nx*3)
      real(kind=8)::xtmp,ytmp,ztmp,r0,tol
      zlo =  0.0
      zhi = 10.56
      r0 = 50
      nzhi = 0
      nzlo = 0
      FEM_zlo = 0.0
      FEM_zhi = 0.0
      tol = 0.01
      do i = 1,nx
        xtmp = g0xyz(i*3-2)
        ytmp = g0xyz(i*3-1)
        ztmp = g0xyz(i*3-0)
        if((abs(ztmp- zhi)<tol).and.
     &     (xtmp*xtmp+ytmp*ytmp)>(r0*r0)) then
           FEM_zhi = FEM_zhi + gpxyz(i*3)-g0xyz(i*3)
           nzhi = nzhi + 1
!           if(me==0) then 
!         write(*,*) 'zhi is found',i,FEM_zhi,nzhi,gpxyz(i*3),g0xyz(i*3)
!           endif
           endif
        if((abs(ztmp-zlo)<tol).and.
     &      (xtmp*xtmp+ytmp*ytmp>r0*r0)) then
!           if(me==0) then 
!              write(*,*) 'zlo is found',i,g0xyz(i*3),nzlo
!           endif
           nzlo = nzlo + 1
           FEM_zlo = FEM_zlo + gpxyz(i*3)-g0xyz(i*3)
         endif
      enddo
      FEM_zhi = FEM_zhi / nzhi
      FEM_zlo = FEM_zhi / nzlo
      if(me==0) then
         write(*,*) 'FEM new size is',FEM_zhi,FEM_zlo,nzhi,nzlo
      endif
      return
      end

      subroutine stiff_media_calc(stress,strain,stiff)
      implicit none
      integer:: i
      real(8):: stress(6),strain(6),stiff(6,6) 
      real(8):: Ls(6)
      real(8):: a(3,2),b(3)
      real(8):: c11,c12,c44
      real(8):: tmp
      real(8):: K(2,2),U(2),Kinv(2,2)
      
      ! use ave_stress, ave_strain to get stiff_media, with independent variables c11 c12 c44
      ! start from c11 + c12 since there is no shear load applied
    
      ! strain(1)*c11 + strain(2)*c21 + strain(3)*c31 = stress(1)
      ! strain(1)*c12 + strain(2)*c22 + strain(3)*c32 = stress(2)
      ! strain(1)*c13 + strain(2)*c23 + strain(3)*c33 = stress(3)
      ! careful, when there is damage, c12 = c21 no longer hold, 
      ! neither does c13 = c31, nor c32 = c23
             

      ! its equation | a11 a12 |               | b1 |            
      !              | a21 a22 | * | x1 x2 | = | b2 | -> A x = b
      !              | a31 a33 |               | b3 |
 
      ! x1 -> c11 , x2 -> c12 = c13, bi -> stress(i)
      ! a11 = strain(1), a12 = strain(2)+strain(3)
      ! a21 = strain(2), a22 = strain(3)+strain(1)
      ! a31 = strain(3), a32 = strain(1)+strain(2)

      ! using least-square method, we got
      !              | a b | * | x1 x2 | = | e |  
      !              | c d |               | f |
      ! where a = a11^2 + a21^2 + a31^1, b = c = a11*a12 + a21*a22 + a31*a32, d = a12^1 + a22^2 + a32^2
      ! and e = a11*b1 + a21*b2 + a31*b3, f = a12*b1 + a22b2 + a32*b3
      
      ! effectively  (A'A) x = (A'b)


      ! therefore 
      do i = 1,6
         Ls(i) = strain(i)*1.0d4
      enddo
      a(1,1) = Ls(1)
      a(2,1) = Ls(2)
      a(3,1) = Ls(3)
      a(1,2) = Ls(2) + Ls(3)
      a(2,2) = Ls(3) + Ls(1)
      a(3,2) = Ls(1) + Ls(2)

      K(1,1) = a(1,1)*a(1,1) + a(2,1)*a(2,1) + a(3,1)*a(3,1)
      K(1,2) = a(1,1)*a(1,2) + a(2,1)*a(2,2) + a(3,1)*a(3,2)
      K(2,1) = a(1,1)*a(1,2) + a(2,1)*a(2,2) + a(3,1)*a(3,2)
      K(2,2) = a(1,2)*a(1,2) + a(2,2)*a(2,2) + a(3,2)*a(3,2)

      U(1) = a(1,1)*b(1)+a(2,1)*b(2)+a(3,1)*b(3)
      U(2) = a(1,2)*b(1)+a(2,2)*b(2)+a(3,2)*b(3)
 
      tmp = K(1,1)*K(2,2) - K(1,2)*K(2,1)    

      if(abs(tmp) < 1e-5) then
         write(*,*) 'error: stiff goes to infinity',tmp
      else
         c11 = ( K(2,2)*U(1) - K(1,2)*U(2))/tmp*1.0d8
         c12 = (-K(2,1)*U(1) + K(1,1)*U(2))/tmp*1.0d8
      endif
      stiff(1,1) = c11
      stiff(2,2) = c11
      stiff(3,3) = c11
      stiff(1,2) = c12
      stiff(1,2) = c12
      stiff(1,3) = c12
      stiff(3,1) = c12
      stiff(2,3) = c12
      stiff(3,2) = c12
      return
      end

      subroutine norm_vector(vector,norm,vector_size)
      implicit none
      integer:: vector_size,i
      real(8):: norm,norm2
      real(8):: vector(vector_size)
      norm = 0.d0
      norm2 = 0.d0
      if(vector_size > 0)then
         do i = 1,vector_size
            norm2 = norm2 + vector(i)*vector(i)
         enddo
         norm = sqrt(norm2)
      endif
      return
      end

