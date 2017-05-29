
!***********************************************************************
!c     Initializes size of modified sparse row format
!***********************************************************************
      subroutine init_size_slu(ijk,ndf,iproc_split,N_update,in_st,isize)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      common/elemnt/nx,nelx,neq,nplane,node
      common/nodeconnec/ielemno,nindx
       dimension ijk(mnelx),ielemno(maxnodelem*maxnode),nindx(maxnode)
       dimension incomp(maxnode)
!***********************************************************************
! c     ijk(I)
! c     N_update(I)
! c     ndf(I)
!***********************************************************************

       indx=1

       do i=1,N_update

          nno=(in_st+i-1)/ndf+1
          icomp=0
          
          if(nno.eq.1)then
            istart=1
          else
            istart=nindx(nno-1)
          endif

          iend=nindx(nno)-1
          do  j=istart,iend

              ielno=ielemno(j)

              do k=1,node

                 imatch=0
                 nn=ijk((ielno-1)*node+k)

                 do ll=1,icomp
                    if(nn.eq.incomp(ll))then
                       imatch=1
                       goto 30
                    endif
                 enddo

  30             if(imatch.eq.0)then

                    icomp=icomp+1
                    incomp(icomp)=nn

                 endif

             enddo
          enddo



          do j=1,nx
             do k=1,icomp
                if(j.eq.incomp(k))then
                   do l=1,ndf
                      indx=indx+1
                   enddo
                endif
             enddo
          enddo


      enddo

      isize=indx-1

      end
!*************************END*******************************************
!***********************************************************************
!c     Initializes modified sparse row format
!***********************************************************************
      subroutine init_msr_slu(ijk,ndf,iproc_split,N_update,in_st,
     & irowptr,ibindx,isize)
!      implicit real*8(a-h,o-z)
      implicit none
      include 'pardis.h'
      

      integer:: ijk(mnelx),ielemno(maxnodelem*maxnode),nindx(maxnode)
      integer:: incomp(maxnode)
!***********************************************************************
! c     ijk(I)
! c     N_update(I)
! c     iupdate(I)-contains rows to be updated by this processor
! c     ibindx(O)-array containing non-zero column indices for each row
! c     consistent with aztec
! c     val(O)-array with elements initialized to zero
! c     ndf(I)
!***********************************************************************
      integer:: ibindx(isize),irowptr(N_update+1)
      integer:: N_update,imatch,iend,ielno
      integer:: icomp,nno,k
      integer:: istart, nn,indx,ll
      integer:: in_st,ndf,i,node
      integer:: j,nx,l,iproc_split
      integer:: isize,nelx,neq,nplane 
      
      
      
      
      common/elemnt/nx,nelx,neq,nplane,node
      common/nodeconnec/ielemno,nindx
      
      indx=1
      
      do i=1,N_update

         nno=(in_st+i-1)/ndf+1
         icomp=0
         irowptr(i)=indx
          
          if(nno.eq.1)then
            istart=1
          else
            istart=nindx(nno-1)
          endif

          iend=nindx(nno)-1
          do  j=istart,iend

              ielno=ielemno(j)

              do k=1,node

                 imatch=0
                 nn=ijk((ielno-1)*node+k)

                 do ll=1,icomp
                    if(nn.eq.incomp(ll))then
                       imatch=1
                       goto 30
                    endif
                 enddo

  30             if(imatch.eq.0)then

                    icomp=icomp+1
                    incomp(icomp)=nn

                 endif

             enddo
          enddo



          do j=1,nx
             do k=1,icomp
                if(j.eq.incomp(k))then
                   do l=1,ndf
                      ibindx(indx)=(j-1)*ndf+l
                      indx=indx+1
                   enddo
                endif
             enddo
          enddo


      enddo
      

      irowptr(N_update+1)=indx-1

      end
      
      
      

!***********************************************************************
!c     Initializes size of modified sparse row format
!***********************************************************************
      subroutine init_size_slu_pbc(ijk,ndf,iproc_split,N_update,in_st,
     & isize,npairs_pbc,id_node_PBC)
      implicit real*8(a-h,o-z)
      include 'pardis.h'
      common/elemnt/nx,nelx,neq,nplane,node
      common/nodeconnec/ielemno,nindx
      dimension ijk(mnelx),ielemno(maxnodelem*maxnode),nindx(maxnode)
      dimension incomp(maxnode)

!***********************************************************************
! c     ijk(I)
! c     N_update(I)
! c     ndf(I)
!***********************************************************************

       indx=1

       do i=1,N_update

          nno=(in_st+i-1)/ndf+1
          icomp=0
          if(nno.eq.1)then
            istart=1
          else
            istart=nindx(nno-1)
          endif

          iend=nindx(nno)-1
          do j=istart,iend
            ielno=ielemno(j)
            do k=1,node
              imatch=0
              nn=ijk((ielno-1)*node+k)
              do ll=1,icomp
                if(nn.eq.incomp(ll))then
                  imatch=1
                  goto 30
                endif
              enddo

  30          if(imatch.eq.0)then
                icomp=icomp+1
                incomp(icomp)=nn
              endif
            enddo
          enddo

         !---------- if PBC ------------!
          kk = 3
          do i_PBC = 1,npairs_PBC(kk)
            id_node1 = id_node_PBC(i_PBC,1,kk)
            id_node2 = id_node_PBC(i_PBC,2,kk)
            if(nno == id_node1)then
              imatch = 0
              do ll = 1,icomp
                if(id_node2 == incomp(ll))then
                  imatch = 1
                  continue
                endif
              enddo
              if(imatch == 0)then
                icomp = icomp + 1
                incomp(icomp) = id_node2
              endif
            elseif(nno == id_node2)then
              imatch = 0
              do ll = 1,icomp
                if(id_node1 == incomp(ll))then
                  imatch = 1
                  continue
                endif
              enddo
              if(imatch == 0)then
                icomp = icomp + 1
                incomp(icomp) = id_node1
              endif
            endif
         enddo
         !---------- PBC end -----------!


          do j=1,nx
             do k=1,icomp
                if(j.eq.incomp(k))then
                   do l=1,ndf
                      indx=indx+1
                   enddo
                endif
             enddo
          enddo


      enddo

      isize=indx-1

      end
!*************************END*******************************************
!***********************************************************************
!c     Initializes modified sparse row format
!***********************************************************************
      subroutine init_msr_slu_pbc(ijk,ndf,iproc_split,N_update,in_st,
     &           irowptr,ibindx,isize,
     &           npairs_pbc,id_node_pbc)
!      implicit real*8(a-h,o-z)
      implicit none
      include 'pardis.h'

      integer:: ijk(mnelx),ielemno(maxnodelem*maxnode),nindx(maxnode)
      integer:: incomp(maxnode)
!***********************************************************************
! c     ijk(I)
! c     N_update(I)
! c     iupdate(I)-contains rows to be updated by this processor
! c     ibindx(O)-array containing non-zero column indices for each row
! c     consistent with aztec
! c     val(O)-array with elements initialized to zero
! c     ndf(I)
!***********************************************************************
      integer:: ibindx(isize),irowptr(N_update+1)
      integer:: N_update,imatch,iend,ielno
      integer:: icomp,nno,k
      integer:: istart, nn,indx,ll
      integer:: in_st,ndf,i,node
      integer:: j,nx,l,iproc_split
      integer:: isize,nelx,neq,nplane 
      integer:: flag_connect
      !---- pbc part -----!
      integer:: i_pbc,ii,jj,kk
      integer:: id_node1,id_node2,id_dof1,id_dof2
      integer:: i_st,i_nd
      common/elemnt/nx,nelx,neq,nplane,node
      common/nodeconnec/ielemno,nindx
      
      indx=1
      
      do i=1,N_update

         nno=(in_st+i-1)/ndf+1
         icomp=0
         irowptr(i)=indx
          
         if(nno.eq.1)then
           istart=1
         else
           istart=nindx(nno-1)
         endif

         iend=nindx(nno)-1
         do j=istart,iend
           ielno=ielemno(j)
           do k=1,node
             imatch=0
             nn=ijk((ielno-1)*node+k)

             do ll=1,icomp
               if(nn.eq.incomp(ll))then
                 imatch=1
                 goto 30
               endif
             enddo

  30         if(imatch.eq.0)then
               icomp=icomp+1
               incomp(icomp)=nn
             endif
           enddo
         enddo

         !---------- if PBC ------------!
         kk = 3
         do i_PBC = 1,npairs_PBC(kk)
            id_node1 = id_node_PBC(i_PBC,1,kk)
            id_node2 = id_node_PBC(i_PBC,2,kk)
            if(nno == id_node1)then
              imatch = 0
              do ll = 1,icomp
                if(id_node2 == incomp(ll))then
                  imatch = 1
                  continue
                endif
              enddo
              if(imatch == 0)then
                icomp = icomp + 1
                incomp(icomp) = id_node2
              endif
            elseif(nno == id_node2)then
              imatch = 0
              do ll = 1,icomp
                if(id_node1 == incomp(ll))then
                  imatch = 1
                  continue
                endif
              enddo
              if(imatch == 0)then
                icomp = icomp + 1
                incomp(icomp) = id_node1
              endif
            endif
         enddo
         !---------- PBC end -----------!

         do j=1,nx
           do k=1,icomp
             if(j.eq.incomp(k))then
               do l=1,ndf
                 ibindx(indx)=(j-1)*ndf+l
                 indx=indx+1
               enddo
             endif
           enddo
        enddo



      enddo
     
      

      irowptr(N_update+1)=indx-1


      end
