      subroutine c_infw( &
                 mot,imot,nmot, &
                 x,y,z,v,mut,tnte1,utau, &
                 vdual,vdual1,vdual2)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action infw.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use sortiefichier
      use modeleturb
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      real mut
!
      dimension imot(nmx)
      dimension ldom(nobj)
      dimension x(ip21),y(ip21),z(ip21)
      dimension v(ip11,ip60)
      dimension mut  (ip12)
      dimension tnte1(ip11,ip60)
      dimension utau(ip42)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60)
!
      call tcmd_infw( &
                 mot,imot,nmot, &
                 ldom,ldomd,kina)
!
      if (kimp.ge.1) then
         call b1_infw(ldom,ldomd,kina)
      endif
!
      do l=1,ldomd
!
        call infw( &
                 ldom(l),x,y,z,v,mut,tnte1, &
                 kina,utau, &
                 vdual,vdual1,vdual2)
!
      enddo
!
      return
      end
