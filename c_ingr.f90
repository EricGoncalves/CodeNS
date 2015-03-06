      subroutine c_ingr( &
                 mot,imot,nmot, &
                 x,y,z)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action ingr.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
	  use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *32 mot(nmx)
!
      dimension imot(nmx)
      dimension ldom(nobj)
      dimension x(ip21),y(ip21),z(ip21)
!
      call tcmd_ingr( &
                 mot,imot,nmot, &
                 ldom,ldomd,king)
!
      if(kimp.ge.1) then
         call b1_ingr(ldom,ldomd,king)
      endif
!
      do l=1,ldomd
        call ingr(ldom(l),x,y,z,king)
      enddo
!
      return
      end
