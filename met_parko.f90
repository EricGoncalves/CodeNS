module mod_met_parko
implicit none
contains
      subroutine met_parko( &
                 l,ncin,ncbd, &
                 v,mut,mu,dist,mnpar)
!
!***********************************************************************
!
!     ACT
!_A   Condition limite a la paroi sur omega pour modeles de Wilcox ou Menter.
!_A   Annulation de la dissipation artificielle pres de la paroi.
!_A   si l=-1: application 1/y^2 a toutes les parois
!_A   si l=0 : application 1/y^2 a toutes les parois
!_A   l>0    : application 1/y^2 uniquement aux parois appartenant au domaine l
!_A            dissipation inchangee
!_A
!_A   LIMITATION : pas de coin. assez de points suivant normale a la paroi
!_A                imxclko superieur ou egal a 1
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    d          : arg real(ip11,ip60 ) ; terme de dissipation artificielle
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    dist
!_I    mnpar
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!_I    nbdko      : com int              ; nombre de parois a traiter
!_I                                        defini dans "atcaldist"
!_I    lbdko      : com int (mtt       ) ; numero interne des parois
!_I                                        (necessaire pour k-omega)
!_I    imxclko    : com int              ; nb pts perpendiculairement paroi
!_I                                        pour cd limite sur omega
!     LOC
!
!***********************************************************************
!
!        isens3=1                             isens3=-1
!
!                  ^                                  ^ //paroi//
!     m3mxf        |---------            m3min        |---------
!       ^          |                                  |  x  nmin
!       |          | . . . .             m3min+isens3 |---------
!       |          |  x                    |      |---|  x  nminp
!       |          |---------              |      |   |---------
!       |      |---|  x  nmaxp             |      |   |  x
!     boucle   |   |---------              |     pond | . . . . .
!       |      |   |  x                    |      |   |  x
!       |     pond | . . . . .           boucle   |   |---------
!       |      |   |  x                    |      |---|  x  nmaxp
!       |      |   |---------              |          |---------
!       |      |---|  x  nminp             |          |  x
!     m3min+isens3 |---------              |          | . . . .
!                  |  x  nmin              v          |
!     m3min        |---------            m3mxf        |---------
!                    //paroi//
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use modeleturb
      use chainecarac
      use sortiefichier
use mod_atindnor
implicit none
integer :: l
integer :: ncin
integer :: ncbd
double precision :: v
double precision :: dist
integer :: mnpar
integer :: isens3
integer :: lpar
integer :: m
integer :: m10
integer :: m1max
integer :: m1min
integer :: m20
integer :: m2max
integer :: m2min
integer :: m3
integer :: m30
integer :: m3clpar
integer :: m3max
integer :: m3min
integer :: m3mxf
integer :: m3tot
integer :: mfbi
integer :: mpar
integer :: mpb0
integer :: mpn0
integer :: mpn1
integer :: n
integer :: ncfic
integer :: nfac
integer :: nmaxp
integer :: nmin
integer :: nminp
integer :: npar
double precision :: pond
double precision :: romegal
!
!-----------------------------------------------------------------------
!
      integer dm1,dm2,dm3
      real mu,mut
      dimension v(ip11,ip60)
      dimension mut(ip12),mu(ip12)
      dimension dist(ip12),mnpar(ip12)
      dimension ncin(ip41),ncbd(ip41)
!
      if(equatt(1:3).eq.'2KO') then
       imxclko=3
       beta1=0.075
      endif
!
!     --------------parametre en dur--------------------------------
!
!     m3clpar : nb mini de points pour traitement condition limite pres paroi
!
      m3clpar=10
!
!     --------------traitement omega en 1/y^2-----------------------
!
      do mpar=1,nbdko
!       boucle sur les parois
!
        mfbi=lbdko(mpar)
        lpar=ndlb(mfbi)
        if(l.le.0 .or. l.eq.lpar) then
!
!         application de la valeur analytique de omega
!
          mpn0=mpn(mfbi)
          mpn1=mpn0+mmb(mfbi)
          mpb0=mpb(mfbi)
!
!                          mpn0 : pointeur dernier point frontiere precedente
!                                 dans tableaux frontieres a normales stockees
!                          mpn1 : pointeur dernier point de la frontiere
!
!                          nmin : pointeur cellule interieure adjacente
!                                 a la paroi
!
!com      calculs des indices direction parallele et normale a la paroi
          call atindnor( &
                 mfbi, &
                 m10,m20,m30, &
                 m1min,m1max,m2min,m2max,m3min,m3max, &
                 dm1,dm2,dm3,isens3)
          m3tot=abs(m3max-m3min)+1
!
          if(m3tot.lt.m3clpar) then
            write(imp,'(/,"!!!met_parko: il faut au moins m3clpar=",i3,' &
            //'"  perpendiculairement a la paroi",/,13x,' &
            //'"domaine =",i3,4x,"nb pts=",i3,4x,"ARRET")')lpar,m3clpar
            stop
          end if
!
!                m3mxf : action du limiteur sur omega entre la paroi et "m3mxf"
!
          m3mxf=m3min+isens3*min(m3clpar,m3tot)
!
          do m=1,mmb(mfbi)
!           boucle sur les points de la frontiere
!
            nfac =mpb0+m
            nmin =ncin(nfac)
!
!                nminp : pointeur second point interieur adjacent a la paroi
!                nmaxp : pointeur dernier point pour ponderation
!
            nminp=nmin+dm3*isens3
            nmaxp=nmin+dm3*isens3*(imxclko+1)
!
!           valeur analytique au point interieur adjacent a la paroi
            romegal=6.*mu(nmin)/(beta1*dist(nmin)**2)
            v(nmin,7)=romegal
!
            do m3=m3min+isens3,m3mxf,isens3
!             boucle sur les "m3mxf" points pres de la paroi
!
!             n=nmin+(m3-m3min)*dm3*isens3
              n=nmin+(m3-m3min)*dm3
              npar=mnpar(n)
              if(npar.gt.mpn0 .and. npar.le.mpn1) then
!
!               la cellule est rattachee a la frontiere
!               A FAIRE : ponderation autorisee si l'ensemble des cellules
!               -------   est rattachee a la frontiere
!
                romegal=6.*mu(n)/(beta1*dist(n)**2)
                if((n-nmaxp)*(n-nminp).le.0) then
                  pond=(dist(n)-dist(nminp))/(dist(nmaxp)-dist(nminp))
                  v(n,7)=romegal*(1.-pond)+pond*v(n,7)
                endif
                v(n,7)=max(v(n,7),romegal)
              endif
!            fin de boucle sur les points pres de la paroi
            enddo
!
!           cellule fictive a la paroi - en principe, toujours rattachee
            ncfic     =ncbd(nfac)
            v(ncfic,7)=v(nmin,7)
!           fin de boucle sur les points de la frontiere
          enddo
        endif
!
      enddo
!
      return
      end
end module
