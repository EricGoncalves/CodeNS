      subroutine svgr( &
                 l,x,y,z, &
                 tn1,tn2,tn3, &
                 disc)
!
!***********************************************************************
!
!     ACT
!_A    Ecriture des coordonnees (x, y, z)
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    disc       : arg char             ; changement de discretisation (centre/noeud)
!_I    kdgv       : com int              ; unite logiq, coordonnees des noeuds
!
!     LOC
!_L    tn1        : arg real(ip00      ) ; tableau de travail
!_L    tn2        : arg real(ip00      ) ; tableau de travail
!_L    tn3        : arg real(ip00      ) ; tableau de travail
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character(len=4 ) :: disc
      dimension x(ip21),y(ip21),z(ip21)
      dimension tn1(ip00),tn2(ip00),tn3(ip00)
!
!     changement eventuel de discretisation
!
      call chdgcv( &
                 l,disc, &
                 x,y,z, &
                 tn1,tn2,tn3, &
                 imin,imax,jmin,jmax,kmin,kmax,kdg)
!
!     stockage du maillage
!
      call writdg( &
                 l,kdg, &
                 imin,imax,jmin,jmax,kmin,kmax, &
                 tn1,tn2,tn3)
!
      return
      end
