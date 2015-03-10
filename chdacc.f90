      subroutine chdacc( &
                 l,disc,ncin,ncbd, &
                 equat, &
                 v,mut, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8, &
                 imin,imax,jmin,jmax,kmin,kmax,kda)
!
!***********************************************************************
!
!     ACT
!_A    Changement eventuel de discretisation pour les variables aerodynamiques
!_A    a partir de valeurs connues aux centres des mailles.
!
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    disc       : arg char             ; changement de discretisation (centre/noeud)
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    equat      : arg char             ; type d'equations modelisant l'ecoulement
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    kdav       : com int              ; unite logiq, var aero aux noeuds
!_I    kdac       : com int              ; unite logiq, var aero aux centres
!_I    kdacf      : com int              ; unite logiq, var aero aux centres
!_I                                        et centres des facettes front
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!
!     OUT
!_O    tn1        : arg real(ip00      ) ; variable de calcul numero 1
!_O    tn2        : arg real(ip00      ) ; variable de calcul numero 2
!_O    tn3        : arg real(ip00      ) ; variable de calcul numero 3
!_O    tn4        : arg real(ip00      ) ; variable de calcul numero 4
!_O    tn5        : arg real(ip00      ) ; variable de calcul numero 5
!_O    tn6        : arg real(ip00      ) ; viscosite turbulente
!_O    imin       : arg int              ; indice min en i
!_O    imax       : arg int              ; indice max en i
!_O    jmin       : arg int              ; indice min en j
!_O    jmax       : arg int              ; indice max en j
!_O    kmin       : arg int              ; indice min en k
!_O    kmax       : arg int              ; indice max en k
!_O    kda        : arg int              ; unite logique, variables
!
!     LOC
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_L    lbd        : com int (mtt       ) ; numero de front a traiter
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use maillage
   use boundary
   use sortiefichier
!
!-----------------------------------------------------------------------
!
      real mut
      character(len=7 ) :: equat
      character(len=4 ) :: disc
      dimension v(ip11,ip60)
      dimension mut(ip12)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00), &
                tn6(ip00),tn7(ip00),tn8(ip00)
      dimension ncbd(ip41),ncin(ip41)
!
      if(disc.eq.'cccc') then
            call cccca( &
                 l,equat, &
                 v,mut, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
      imin=ii1(l)
      imax=ii2(l)-1
      jmin=jj1(l)
      jmax=jj2(l)-1
      kmin=kk1(l)
      kmax=kk2(l)-1
      kda=kdac
      elseif(disc.eq.'cccf') then
            call cccca( &
                 l,equat, &
                 v,mut, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
      imin=ii1(l)-1
      imax=ii2(l)
      jmin=jj1(l)-1
      jmax=jj2(l)
      kmin=kk1(l)-1
      kmax=kk2(l)
      kda=kdacf
!
      elseif(disc.eq.'cccv') then
       nbd=0
       do mfb=1,mtbx
        if(ndlb(mfb).eq.l) then
         nbd=nbd+1
         lbd(nbd)=mfb
        endif
       enddo
       call rfspstf(v(1,1),ncin,ncbd)
       call rfspstc(l,v(1,1))
       call rfspstf(v(1,2),ncin,ncbd)
       call rfspstc(l,v(1,2))
       call rfspstf(v(1,3),ncin,ncbd)
       call rfspstc(l,v(1,3))
       call rfspstf(v(1,4),ncin,ncbd)
       call rfspstc(l,v(1,4))
       call rfspstf(v(1,5),ncin,ncbd)
       call rfspstc(l,v(1,5))
       if(equat(1:2).eq.'ns') then
         call rfspstf(mut,ncin,ncbd)
         call rfspstc(l,mut)
       endif
       if(equat(6:7).eq.'ke') then
         call rfspstf(v(1,6),ncin,ncbd)
         call rfspstc(l,v(1,6))
         call rfspstf(v(1,7),ncin,ncbd)
         call rfspstc(l,v(1,7))
       endif
        call cccva( &
                 l,equat, &
                 v,mut, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
        imin=ii1(l)
        imax=ii2(l)
        jmin=jj1(l)
        jmax=jj2(l)
        kmin=kk1(l)
        kmax=kk2(l)
        kda=kdav
      endif
!
      return
      end
