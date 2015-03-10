module mod_chrono
implicit none
contains
      subroutine chrono( &
                 img, &
                 mu,mut,u,dt,dtmin, &
                 sn,vol, &
                 tn1,tn2,tn3,tn4, &
                 cson,pression)
!
!***********************************************************************
!
!_DA  DATE_C : octobre 2006 - Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul du pas de temps pour tous les points de la configuration
!_A    de calcul.
!
!
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    t          : arg real(ip11,ip60 ) ; variables de calcul
!_I    equat      : arg char             ; type d'equations modelisant l'ecoule ment
!_I    sn         : arg real(lgsnlt,nind,ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    lgsnlt     : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    gam1       : com real             ; rap chal spec -1
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
!
!     OUT
!_O    dt         : arg real(ip11      ) ; pas de temps
!
!     LOC
!_L    sfsi       : arg real(ip00      ) ; surface au carre d'une facette i
!_L    sfsj       : arg real(ip00      ) ; surface au carre d'une facette j
!_L    sfsk       : arg real(ip00      ) ; surface au carre d'une facette k
!_L    dism       : arg real(ip00      ) ; distance caracteristique d'une maille
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use schemanum
      use constantes
      use chainecarac
      use sortiefichier 
use mod_chronos_prcd
use mod_chronos
implicit none
integer :: indc
integer :: img
double precision :: u
double precision :: dt
double precision :: sn
double precision :: vol
double precision :: tn1
double precision :: tn2
double precision :: tn3
double precision :: tn4
double precision :: cson
double precision :: pression
integer :: i
integer :: j
integer :: k
integer :: l
integer :: lgsnlt
integer :: lm
integer :: n
integer :: ndeb
integer :: nfin
integer :: nid
integer :: nijd
integer :: njd
integer :: npsn
!
!-----------------------------------------------------------------------
!
      real mu,mut,dtmin
      dimension u(ip11,ip60)
      dimension dt(ip11),cson(ip11),vol(ip11),pression(ip11)
      dimension mu(ip12),mut(ip12)
      dimension sn(ip31*ndir)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00)

      indc(i,j,k)=npc(lm)+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd

!-----calcul du pas de temps local en stationnaire--------------------
!
      do l=1,lzx
       lm=l+(img-1)*lz
       npsn  =ndir*npfb(lm)+1
       lgsnlt=nnn(lm)
!
       if(kprec.eq.0) then
            call chronos( &
                 lm,eta(lm),mu,mut, &
                 u, &
                 dt,equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 tn1,tn2,tn3,tn4, &
                 cson)
!
      elseif(kprec.ge.1) then
           call chronos_prcd( &
                 lm,eta(lm),mu,mut, &
                 u, &
                 dt,equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 tn1,tn2,tn3,tn4, &
                 cson,pression)
      endif
!
      enddo
!
!-----calcul instationnaire avec pas de temps global---------------
!
      if(kdtl.eq.0) then
       dtmin=reelmx
!
       do l=1,lzx
        lm=l+(img-1)*lz
        nid  = id2(lm)-id1(lm)+1
        njd  = jd2(lm)-jd1(lm)+1
        nijd = nid*njd
        do k = kk1(lm),kk2(lm)-1
         do j = jj1(lm),jj2(lm)-1
          do i = ii1(lm),ii2(lm)-1
           n = indc(i,j,k)
           dtmin=min(dtmin,dt(n))
          enddo
         enddo
        enddo
!
      enddo
!
      dtmin=min(dt1min,dtmin)
      if(dtmin.lt.dt1min) then
       write(imp,'("===> chrono")')
       write(imp,'(2(1pe11.3))') &
           dtmin,dt1min
       STOP
      endif
!
      do l=1,lzx
       lm=l+(img-1)*lz
       ndeb = indc(id1(lm),jd1(lm),kd1(lm))
       nfin = indc(id2(lm),jd2(lm),kd2(lm))
       do n = ndeb,nfin
        dt(n)=dtmin
       enddo
      enddo
!
      endif
!
      return
      end
end module
