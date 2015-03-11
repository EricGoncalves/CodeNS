module mod_atcaldis
implicit none
contains
      subroutine atcaldis( &
                 x,y,z,nxn,nyn,nzn, &
                 xpar,ypar,zpar,xcc,ycc,zcc,dist2, &
                 dist,mnpar,fgam,img, &
                 ncin,mnc,ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Calcul de la distance caracteristique de chaque cellule a une
!_A    paroi et rattachement de chaque cellule a une facette de paroi.
!_A
!      REMARQUE: les tableaux xpar, ypar, zpar, dist2 sont dimensionnes
!                sur le plus grand domaine (ip00) alors qu'ils devraient
!                l'etre sur les facettes frontieres a normales stockees (ip42).
!                On verifie que ip00 est superieur a ip42.
!
!_I    ncycle     : arg int              ; nbr tot de cycles de l'execution
!_I                                        courante
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    exs1       : arg real             ; premier coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    exs2       : arg real             ; deuxieme coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    reelmx     : com real             ; nombre reel grand
!_I    lzx        : com int              ; nbr total de domaines
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mtnx       : com int              ; nbr total de frontieres  a normales stockes
!_I    mtcx       : com int              ; nbr total de frontieres coincidentes
!_I    mtrx       : com int              ; nbr total de frontieres recouvertes
!_I    mtax       : com int              ; nbr total de frontieres autres
!_L    ndimctbx   : com int              ; nbr de cellules de tts les domaines
!_L                                        (pts fictifs inclus)
!_L    mdimtnx    : com int              ; nbr de pts de ttes les front
!_L                                        a normales stockees
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    nnc        : com int (lt        ) ; nombre de cellules du dom (dont fic.)
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
!_I    kd2        : com int (lt        ) ; indice max en k fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_L    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_L    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_L                                        en fct du numero externe
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    nfbr       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere recouverte
!_I    nfbc       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere coincidente
!_I    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une front a normales stockees
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    indfl      : arg char             ; type de plan de la frontiere
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    kcaldis    : com int              ; cle calcul distance aux parois
!_I    reelmx     : com real             ; nombre reel grand
!
!     OUT
!_O    dist       : arg real(ip12      ) ; distance a la paroi
!_O    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_O                                        stockees du point de rattach normale
!_O    nbdko      : com int              ; nombre de parois a traiter
!_O    lbdko      : com int (mtt       ) ; numero interne des parois
!_O                                        (necessaire pour k-omega)
!
!     LOC
!_L    tn1        : arg real(ip00      ) ; tableau de travail -> xpar
!_L    tn2        : arg real(ip00      ) ; tableau de travail -> ypar
!_L    tn3        : arg real(ip00      ) ; tableau de travail -> zpar
!_L    tn4        : arg real(ip00      ) ; tableau de travail -> xcc
!_L    tn5        : arg real(ip00      ) ; tableau de travail -> ycc
!_L    tn6        : arg real(ip00      ) ; tableau de travail -> zcc
!_L    tn7        : arg real(ip00      ) ; tableau de travail ->dist2
!_L    xpar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_L    ypar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_L    zpar       : arg real(ip00      ) ; coordonnee y centre facette paroi
!_L    xcc        : arg real(ip00      ) ; coordonnee z centre cellule
!_L    ycc        : arg real(ip00      ) ; coordonnee y centre cellule
!_L    zcc        : arg real(ip00      ) ; coordonnee z centre cellule
!_L    dist2      : arg real(ip00      ) ; travail
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_L    lbd        : com int (mtt       ) ; numero interne de front a traiter
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use sortiefichier
   use maillage
   use constantes
   use boundary
use mod_atdist_2
use mod_atccfp
use mod_atdist_1
use mod_at_ecrdist
use mod_atecrfp
use mod_atccc
use mod_atindnor
implicit none
double precision :: x
double precision :: y
double precision :: z
double precision :: xpar
double precision :: ypar
double precision :: zpar
double precision :: xcc
double precision :: ycc
double precision :: zcc
double precision :: dist2
double precision :: dist
integer :: mnpar
double precision :: fgam
integer :: img
integer :: ncin
integer :: mnc
integer :: ncbd
integer :: ierr
integer :: isens3
integer :: l
integer :: lm
integer :: m
integer :: m10
integer :: m1max
integer :: m1min
integer :: m20
integer :: m2max
integer :: m2min
integer :: m30
integer :: m3max
integer :: m3min
integer :: mf
integer :: n
integer :: nfbe
integer :: nfbi
integer :: no
!
!-----------------------------------------------------------------------
!
      integer dm1,dm2,dm3
      real nxn,nyn,nzn
      dimension dist(ip12),mnpar(ip12),fgam(ip42)
      dimension x(ip21),y(ip21),z(ip21)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension xpar(ip00),ypar(ip00),zpar(ip00)
      dimension xcc (ip00),ycc (ip00),zcc (ip00),dist2(ip00)
      dimension ncin(ip41)
      dimension mnc(ip43)
!
!     -----------------------------------------------------------------
!
!     verification compatibilite du dimensionnement
!     Les centres des facettes parois a normales stockees sont mis
!     dans xpar(ip00), ypar(ip00),zpar(ip00),dist2(ip00). Il faut donc :
!     ndimub=ip00 > mdimtnx=ip42
      if(ip00 .le. mdimtnx) then
        write(imp,'(/,''!!!atcaldis: ip00='',i6,4x,''inferieur a mdimtnx='',i6)')ip00,mdimtnx
        stop
      endif
!
!     -----------------------------------------------------------------
!     Liste des parois - numero interne
!                        nba : rang de traitement
!                        nbd : nombre de frontieres a traiter
!
      ierr =0
      nbd  =0
      nbdko=0
      do no=1,mtbx
        nfbe=nba(no)
        nfbi=nfei(nfbe)
!      modification du test pour la condition aux limites lois de paroi et calcul Euler
        if((cl(nfbi)(1:2).eq.'pa').or.(cl(nfbi)(1:2).eq.'lp').or. &
           (cl(nfbi)(1:2).eq.'gl')) then  !la frontiere est une paroi
!
          nbd         =nbd+1
          lbd(nbd)    =nfbi
          nbdko       =nbdko+1
          lbdko(nbdko)=nfbi
!
!         verification qu'elle est definie a normale stockee
          do m=1,mtnx
             if(nfbn(m).eq.nfbi) then
                ierr=-1
                exit
             endif
          enddo
          ierr=ierr+1
        endif
      enddo
      if(ierr.ne.0) then
        write(imp,'(/,"!!!atcaldis: pb declaration parois - arret")')
        stop
      end if
      if(nbd.eq.0) then
        write(imp,'(/,"!!!atcaldis: pas de paroi pour calcul des distances aux parois - arret")')
        stop
      end if
!
!     ----------------------------------------------------------
!       initialisation des distances a une valeur infinie
!
        do n=1,ndimctbx
          dist(n) =reelmx
          mnpar(n)=0
        enddo
!     ----------------------------------------------------------
!
!     boucle sur les frontieres a traiter (parois)
!
      do mf=1,nbd
        nfbi=lbd(mf)
!
!       calcul des indices pour deplacement "parallelement" a la paroi
        call atindnor( &
                 nfbi, &
                 m10,m20,m30, &
                 m1min,m1max,m2min,m2max,m3min,m3max, &
                 dm1,dm2,dm3,isens3)
!
!       calcul des centres des facettes formant les parois
        call atccfp( &
                 x,y,z, &
                 xpar,ypar,zpar, &
                 ncin, &
                 nfbi,dm1,dm2,dm3,isens3)
!
!       ecriture centre facettes parois
!
        call atecrfp( &
                 x,y,z, &
                 xpar,ypar,zpar, &
                 nxn,nyn,nzn, &
                 ncin, &
                 nfbi)
!
        if(kcaldis.eq.1) then
!
!         integration de la distance suivant les lignes du maillage
!         qui partent de la paroi
!
          l=ndlb(nfbi)
!
!         calcul des centres des cellules
          call atccc( &
                 x,y,z, &
                 xcc,ycc,zcc, &
                 l)
!
!         integration le long du maillage
          call atdist_1( &
                 x,y,z, &
                 nxn,nyn,nzn, &
                 xpar,ypar,zpar, &
                 xcc,ycc,zcc, &
                 dist,mnpar, &
                 ncin,m3min,m3max,dm3,isens3, &
                 nfbi)
        end if
      enddo
!
!     -------------------------------------------------------------
!
      if(kcaldis.eq.2) then

!       Calcul de la distance de toutes les cellules a la facette
!       la plus proche. Toutes les facettes de toutes les parois
!       sont esssayees.
!       boucle sur les domaines
!
        do l=1,lzx
         lm=l+(img-1)*lz
!
!         calcul des centres des cellules
          call atccc( &
                 x,y,z, &
                 xcc,ycc,zcc, &
                 lm)
!
!         calcul distrance de chaque cellule d'un domaine aux parois et
!         cellule origine de la "normale"
          call atdist_2( &
                 x,y,z, &
                 xpar,ypar,zpar, &
                 xcc,ycc,zcc,dist2, &
                 dist,mnpar, &
                 lm)
          if(kecrdis.eq.1) then
!           ecriture disque des distances (fichiers separes "fdist_l")
            call at_ecrdist( &
                 lm,         &
                 dist,mnpar)
          endif
        enddo
      endif
!
      return
      end subroutine
end module
