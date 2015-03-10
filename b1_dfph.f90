      subroutine b1_dfph
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfph'.
!
!     INP
!_I    equat      : com char             ; type d'equations modelisant l'ecoulement
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    rd         : com real             ; constante des gaz parfaits
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
!_I    reynz      : com real             ; nombre de Reynolds calcule avec
!_I                                        les grandeurs d'adimensionnement,
!
!***********************************************************************
!
      use sortiefichier
      use chainecarac
      use kcle
      use proprieteflu
implicit none
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=24) ::  cgam,crd,cpinfl,cql,cpr,cprt,creynz
!
      call convich(kgam,cgam)
      call convich(krd,crd)
      call convich(kpinfl,cpinfl)
      call convich(kql,cql)
      call convich(kpr,cpr)
      call convich(kprt,cprt)
      call convich(kreynz,creynz)
!
       form='(/,2x,''definition de la physique''/' &
             //'2x,''-------------------------'',/' &
             //'2x,''gamma liquide            : '',4x,e12.6,2x,a/' &
             //'2x,''cte gaz liquide          : '',4x,e12.6,2x,a/' &
             //'2x,''pinf liquide             : '',4x,e12.6,2x,a/' &
             //'2x,''energie ref liquide      : '',4x,e12.6,2x,a)'
       write(imp,form) gam,cgam, &
                       rd,crd, &
                       pinfl,cpinfl, &
                       ql,cql
      if (equat(1:2).eq.'ns') then
       form  ='(2x,''pr liquide               : '',4x,e12.6,2x,a/' &
             //'2x,''prt                      : '',4x,e12.6,2x,a/' &
             //'2x,''reynz                    : '',4x,e12.6,2x,a)'
       write(imp,form) pr,cpr, &
                       prt,cprt, &
                       reynz,creynz
      endif
!
      return
      end
