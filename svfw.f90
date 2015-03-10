      subroutine svfw( &
                 l,v,mut,utau, &
                 disc,ncin,ncbd, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
!
!***********************************************************************
!
!     ACT
!_A    Ecriture des variables aerodynamiques (v) (et eventuellement de la
!_A    viscosite turbulente (mut)) selon le type de discretisation specifie.
!
!     INP
!_I    l          : arg int             ; numero de domaine
!_I    v          : arg real(ip11,ip60) ; variables a l'instant n+alpha
!_I    mut        : arg real(ip12     ) ; viscosite turbulente
!_I    disc       : arg char            ; changement de discretisation (centre/noeud)
!_I    ncin       : arg int (ip41     ) ; ind dans un tab tous domaines de la
!_I                                       cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41     ) ; ind dans un tab tous domaines d'une
!_I                                       cellule frontiere fictive
!_I    equat      : com char            ; type d'equations modelisant l'ecoulement
!_I    kdac       : com int             ; unite logiq, var aero aux centres
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use sortiefichier
      use chainecarac
!
!-----------------------------------------------------------------------
!
      character(len=4 ) :: disc
      real mut
!
      dimension utau(ip42)
      dimension v(ip11,ip60)
      dimension mut(ip12)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00), &
                tn5(ip00),tn6(ip00),tn7(ip00),tn8(ip00)
      dimension ncin(ip41),ncbd(ip41)
!
!     changement eventuel de discretisation
!
       call chdacc( &
                 l,disc,ncin,ncbd, &
                 equat, &
                 v,mut, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8, &
                 imin,imax,jmin,jmax,kmin,kmax,kdac)
!
!     stockage des variables aerodynamiques
!
       call writda( &
                 l,kdac,equat,utau, &
                 imin,imax,jmin,jmax,kmin,kmax, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,ip12,tn8)
!
      return
      end
