      subroutine c_inbdn( &
                 mot,imot,nmot, &
                 x,y,z, &
                 sn, &
                 ncbd,nxn,nyn,nzn, &
                 tn1,tn2,tn3,tn4,tn5,tn6, &
                 tn7,tn8,tn9)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action inbdn.
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
      character(len=32) ::  mot(nmx)
!
      real nxn,nyn,nzn
!
      dimension imot(nmx)
      dimension lmfb(mtb)
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension sn(ip31*ndir)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
!
      dimension tn1(ip00),tn2(ip00),tn3(ip00), &
                tn4(ip00),tn5(ip00),tn6(ip00), &
                tn7(ip00),tn8(ip00),tn9(ip00)
!
      call tcmd_inbdn( &
                 mot,imot,nmot, &
                 lmfb,lmfbd,kibdn)
!
      if (kimp.ge.1) then
            call b1_inbdn(lmfb,lmfbd,kibdn)
      endif
!
      do l=1,lmfbd
!
        call inbdn( &
                 lmfb(l),kibdn, &
                 x,y,z, &
                 sn, &
                 ncbd,nxn,nyn,nzn, &
                 tn1,tn2,tn3,tn4,tn5,tn6, &
                 tn7,tn8,tn9)
!
      enddo
!
      return
      end
