      subroutine infw( &
                 l,x,y,z,v,mut,tnte1, &
                 kina,utau, &
                 vdual,vdual1,vdual2)
!
!***********************************************************************
!
!_DA  DATE_C : Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Initialisation des variables aerodynamiques, c'est-a-dire
!_A    des tableaux v et mut: ro ,ro*u ,ro*v ,ro*w, ro*e, mut.
!
!     VAL
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    kina       : arg int              ; cle initialisation var de calcul
!_I    equat      : com char             ; type d'equations modelisant l'ecoulement
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    omg        : com real             ; vitesse rotation du repere relatif
!_I    kdac       : com int              ; unite logiq, var aero aux centres
!
!     OUT
!_O    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use chainecarac
      use definition
      use sortiefichier
      use schemanum
implicit none
integer :: l
double precision :: x
double precision :: y
double precision :: z
double precision :: v
double precision :: tnte1
integer :: kina
double precision :: utau
double precision :: vdual
double precision :: vdual1
double precision :: vdual2
integer :: img
integer :: keinit
integer :: lm
!
!-----------------------------------------------------------------------
!
      real mut
      dimension x(ip21),y(ip21),z(ip21)
      dimension v(ip11,ip60),tnte1(ip11,ip60)
      dimension mut(ip12)
      dimension utau(ip42)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60)
!
      do img=1,lgx
       lm=l+(img-1)*lz
       if (kina.eq.1) then
!      remplissage des tableaux v et mut par un sous-programme d' initialisation
 
        if ((img.eq.1).or.(kfmg.ge.1)) then
          call utinia( &
                 lm,x,y,z,v,mut, &
                 kina,                 &
                 vdual,vdual1,vdual2)
          endif

       elseif(kina.eq.0) then
!      remplissage des tableaux v et mut par lecture d'un fichier de reprise
         if(img.eq.1) then
           call readda( &
                 l,kdac, &
                 v,mut,utau, &
                 vdual,vdual1,vdual2)
!
         elseif(img.gt.1.and.kfmg.ge.1.and.kfmg.lt.3) then
            call utinia( &
                 lm,x,y,z,v,mut, &
                 kina, &
                 vdual,vdual1,vdual2)
         endif

      elseif(kina.eq.2) then
       if(img.eq.1) then
            call utreadav( &
                 l,kdav,equat, &
                 tnte1,v,mut,keinit)
!
        elseif((img.gt.1).and.(kfmg.ge.1)) then
            call utinia( &
                 lm,x,y,z,v,mut, &
                 kina, &
                 vdual,vdual1,vdual2)
        endif
       endif
      enddo
!
      return
      end
