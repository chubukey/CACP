      subroutine   PostProcessing(nx,nelx,nprocs_cpfem,nstep,
     & gpxyz, ijk)
      implicit none
      include 'pardis.h'
      integer:: nx,nelx,nprocs_cpfem,nstep, nstp, dim2
      integer:: elemstr, elemend, dummy, ijk(mnelx)
      integer:: i,j,k,ii,jj,kk,n,t
      real(8):: gpxyz(maxcrd)
      real(8):: displacement(nx,3),fpout(nelx,6)
      integer:: nodeset(4), nodelem(maxnodelem,nx)
      real(8):: nodalvalue(nx,6), totvolume, weight,distance
      real(8):: center(nelx,3), xset(4),yset(4),zset(4)
      real(8):: xpos,ypos,zpos
      character*10 fname
      character (LEN=20) ::  filname
      
      dim2=6
      
      
      
      write(*,*)('****************************************')
      write(*,*)('****************************************')
      write(*,*)('Now start to write tecplot output file')
      
      
      open(601,file='disp.out')
      
      
      nstp=nstep/npost
      
      do t=1,nstp
      do i=1,nx
         read(601,*)(displacement(i,j), j=1,3)           !only last step
      enddo
      enddo
      close(601)
      
      write(*,*)('finish reading displacement')
      
      
      do ii=0,nprocs_cpfem-1
      call get_subscript(ii,fname)      
      filname='stressfs'//fname
      open(602,file=filname)
         do t=1,nstp
            read(602,*)elemstr, elemend
            do i=elemstr,elemend
               do j=1,dim2 
                    read(602,*)fpout(i,j), dummy
                    if (dummy /=i)write(*,*)('error code 1')     !still only record last step
               enddo
            enddo
          enddo  
         close(602)
       enddo
       
      write(*,*)('finish reading stress')
      
      
      
! find for each node, which element it belongs to 
      do n=1,nx
      do i=1,maxnodelem
      nodelem(i,n)=0
      enddo
      enddo
      
      do n=1,nx
       k=0
       kk=0
       do i=1,nelx
            do j=1,4
            kk=kk+1
              if(ijk(kk).eq.n)then
              k=k+1
              nodelem(k,n)=i
              endif
            enddo
        enddo
        
        
        if (k.eq.0)then
        write(*,*)'this node is not connected to any element', n
        endif
        
        if (k.gt.maxnodelem) then
        write(*,*)'maxnodelem is not large enough', k
        endif
        
        enddo



       
! calculte the position of each elem's gauss point
      do k=1,nelx
           nodeset(1)=ijk(k*4-3)
           nodeset(2)=ijk(k*4-2)
           nodeset(3)=ijk(k*4-1)
           nodeset(4)=ijk(k*4)
            do j=1,4
            xset(j)=gpxyz(nodeset(j)*3-2)
            yset(j)=gpxyz(nodeset(j)*3-1)
            zset(j)=gpxyz(nodeset(j)*3)
            enddo

       center(k,1)=1.0/4.0*(xset(1)+xset(2)+xset(3)+xset(4))
       center(k,2)=1.0/4.0*(yset(1)+yset(2)+yset(3)+yset(4))
       center(k,3)=1.0/4.0*(zset(1)+zset(2)+zset(3)+zset(4)) 
          
       enddo    
      
      
       write(*,*)'finish calculating the center of each elem'
      
      
! now calculate nodal value of stress  
        do n=1,nx
        xpos=gpxyz(n*3-2)
        ypos=gpxyz(n*3-1)
        zpos=gpxyz(n*3)
        
        
        do j=1,dim2
        nodalvalue(n,j)=0.0d0
        totvolume=0.0
        do i=1,maxnodelem
          k=nodelem(i,n)
          if (k /= 0) then
         distance=sqrt((xpos-center(k,1))**2+(ypos-center(k,2))**2
     &   +(zpos-center(k,3))**2)
         weight=exp(-3*distance)
         nodalvalue(n,j)=nodalvalue(n,j) + weight*fpout(k,j)
         totvolume=totvolume+weight
          endif
        enddo
         nodalvalue(n,j)=nodalvalue(n,j)/totvolume
         enddo 
      enddo  

        write(*,*) ('finish calculating nodal value') 
      
      
       ! at last, output to tecplot format
      
      open(603, file='TecplotInput_stress.dat')
        
      write(603,*)('VARIABLES=   "X","Y","Z","stress11","stress22",
     & "stress33","stress12","stress13","stress23" ')

      write(603,*) 'ZONE N=  ', nx ,'  E=  ',nelx, 
     &  'F=FEPOINT, ET=TETRAHEDRON'

      do i=1,nx
      write(603,*) gpxyz(i*3-2),gpxyz(i*3-1),gpxyz(i*3),
     & (nodalvalue(i,j),j=1,dim2)
      enddo
      
      do i=1,nelx
      write(603,*) ijk(i*4-3),ijk(i*4-2),ijk(i*4-1),ijk(i*4)
      enddo

      close(603)
      
      write(*,*) ('finish writing output file') 
      
      
      return
      end




