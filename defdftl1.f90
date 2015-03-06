      subroutine defdftl1
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dftl1.
!
!***********************************************************************
!
      use kcle
      use chainecarac
!
!-----------------------------------------------------------------------
!
     do itit=1,80
      titrt1(itit:itit)=' '
     enddo
     titrt1='  CALCUL RANS  '
     ktitrt1=1
!
      return
      end
