module mod_cllparoi1
implicit none
contains
      subroutine cllparoi1( &
                 img,ncyc, &
                 v,dist, &
                 mu,mut, &
                 nxn,nyn,nzn, &
                 ncin,ncbd, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 qcx,qcy,qcz, &
                 mnpar,fgam,tp, &
                 ztemp)
!
!***********************************************************************
!
!_DA  DATE_C :decembre 1998-- AUTEUR : Eric Goncalves
!
!     ACT
!_A    Application de la condition aux limites lois de paroi analytiques.
!
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_I                                        contraintes visqueuses
!_I    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_I                                        contraintes visqueuses
!_I    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_I                                        contraintes visqueuses
!_I    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_I                                        contraintes visqueuses
!_I    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_I                                        contraintes visqueuses
!_I    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_I                                        contraintes visqueuses
!_I    qtx        : arg real(ip40      ) ; composante en x flux de chaleur
!_I    qty        : arg real(ip40      ) ; composante en y flux de chaleur
!_I    qtz        : arg real(ip40      ) ; composante en z flux de chaleur
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    img        : com int              ; niveau de grille (multigrille)
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npfb       : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes facettes
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpd        : com int (mtt       ) ; tableau de pointeur sur l'ensemble
!_I                                        des frontieres
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
use mod_lparoi1
use mod_lparoi2
implicit none
integer :: img
integer :: ncyc
double precision :: v
double precision :: dist
integer :: ncin
integer :: ncbd
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: qcx
double precision :: qcy
double precision :: qcz
integer :: mnpar
double precision :: fgam
double precision :: tp
double precision :: ztemp
integer :: l
integer :: lm
integer :: mfb
integer :: mfbm
integer :: no
!
!-----------------------------------------------------------------------
!
      double precision mu,mut,nxn,nyn,nzn
!
      dimension mu(ip12),mut(ip12), &
                toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12), &
                qcx(ip12),qcy(ip12),qcz(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension ncbd(ip41),ncin(ip41),dist(ip12)
      dimension v(ip11,ip60)
      dimension mnpar(ip12),fgam(ip42)
      dimension tp(ip40)
      dimension ztemp(ip11)
!
!    boucle sur toutes les frontieres
      do no=1,mtbx
!
      mfb=nba(no)
      mfbm=mfb+(img-1)*mtb
      l=ndlb(mfb)
      lm=l+(img-1)*lz
!
!      lois de paroi en parois adiabatiques
!
        if (cl(mfb)(1:3).eq.'lp2') then
          call lparoi1( &
                 lm,ncyc, &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfbm, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 qcx,qcy,qcz, &
                 mnpar,fgam, &
                 ztemp)
!
!      lois de paroi en parois isothermes
!
        else if (cl(mfb)(1:3).eq.'lp3') then
          call lparoi2( &
                 lm,ncyc, &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfbm, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 qcx,qcy,qcz, &
                 mnpar,fgam,tp, &
                 ztemp)
!
        endif
!
      enddo
!
      return
      end subroutine
end module
