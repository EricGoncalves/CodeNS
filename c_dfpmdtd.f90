      subroutine c_dfpmdtd(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfpmdtd.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    pctvort    : com real(lt        ) ; pourcentage de tourbillon pour
!_I                                        calcul d'epaisseur de couche lim
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
      dimension ldom(nobj)
      dimension lgr(nobj)
!
      call tcmd_dfpmdtd( &
                 mot,imot,nmot, &
                 ldom,ldomd, &
                 lgr,lgrd)
!
      if(kimp.ge.1) then
            call b1_dfpmdtd( &
                 ldom,ldomd, &
                 lgr,lgrd)
      endif
!
      return
      end
