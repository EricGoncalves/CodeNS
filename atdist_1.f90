          subroutine atdist_1( &
                 x,y,z, &
                 nxn,nyn,nzn, &
                 xpar,ypar,zpar, &
                 xcc,ycc,zcc, &
                 dist,mnpar, &
                 ncin,m3min,m3max,dm3,isens3, &
                 nfbi)
!
!***********************************************************************
!
!     ACT
!_A    calcul de la distance du centre d'une cellule a la paroi
!_A    en integrant suivant la ligne de maillage issue de la paroi.
!_A    La distance est croissante et projetee sur la normale a la paroi
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    xpar       : arg real(ip00      ) ; coordonnee x centre facette paroi
!_I    ypar       : arg real(ip00      ) ; coordonnee y centre facette paroi
!_I    zpar       : arg real(ip00      ) ; coordonnee z centre facette paroi
!_I    xcc        : arg real(ip00      ) ; coordonnee x centre cellule
!_I    ycc        : arg real(ip00      ) ; coordonnee y centre cellule
!_I    zcc        : arg real(ip00      ) ; coordonnee z centre cellule
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    m3min      : arg int              ; indice debut normale paroi
!_I    m3max      : arg int              ; indice fin normale paroi
!_I    isens3     : arg int              ; sens variation m3 en s'eloignant de
!_I                                        la paroi
!_I    nfbi       : arg int              ; numero interne de la frontiere
!_L    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!
!     OUT
!_O    dist       : arg real(ip12      ) ; distance a la paroi
!_O    mnpar      : arg real(ip12      ) ; pointeur tablaux front normales
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use sortiefichier
   use maillage
   use boundary
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
double precision :: dist
integer :: mnpar
integer :: ncin
integer :: m3min
integer :: m3max
integer :: isens3
integer :: nfbi
double precision :: distn
integer :: ierr
integer :: l
integer :: m
integer :: m0b
integer :: m0n
integer :: m3
integer :: mb
integer :: mbn
integer :: mt
integer :: n0
integer :: nc
integer :: nci
integer :: ncil
integer :: ncl
integer :: ncm1
integer :: ncm1l
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
      integer dm3
      dimension x(ip21),y(ip21),z(ip21)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension xpar(ip00),ypar(ip00),zpar(ip00)
      dimension xcc (ip00),ycc (ip00),zcc (ip00)
      dimension dist(ip12),mnpar(ip12)
      dimension ncin(ip41)
!
      m0b=mpb(nfbi)
      m0n=mpn(nfbi)
      mt =mmb(nfbi)
      l  =ndlb(nfbi)
      n0 =npc(l)
!
!                   m0b mb  : pointeur dans tab toutes frontieres
!                   m0n mbn : pointeur dans tab frontieres a normales stockees
!                   mt      : nombre de points de la frontiere
!                   nci     : pointeur cellule adjacente interieure a la
!                             frontiere, numerotation tous domaines
!                   ncil    : pointeur cellule adjacente interieure a la
!                             frontiere, numerotation locale au domaine
!                   n0      : deniere cellule domaine precedent dans
!                             tableaux toutes cellules tous domaines
!                   nc      : pointeur tous domaines
!                   ncl     : pointeur local au domaine
!                   dist    : distance, tableau tous domaines
!                   xcc     : centre celulle tableau 1 domaine
!
!
!     boucle sur tous les points de la paroi
      do m=1,mt
        ierr=0
        mb =m0b+m
        mbn=m0n+m
!
!       initialisation de l'integration des distances
!       projection de la disatnce sur la normale
!
        nci  =ncin(mb)
        ncil =nci-n0
        distn=(xcc(ncil)-xpar(mbn))*nxn(mbn)+ &
              (ycc(ncil)-ypar(mbn))*nyn(mbn)+ &
              (zcc(ncil)-zpar(mbn))*nzn(mbn)
        if(distn.lt.0.) ierr=ierr+1
        if(distn.lt.dist(nci)) then
          dist(nci) =distn
          mnpar(nci)=mbn
        else !passage au point suivant de la frontiere
          exit
        endif
!
!       integration suivant la ligne de maillage
        ncm1l=ncil
        do m3=m3min+isens3,m3max-isens3,isens3
          ncl  =ncm1l+dm3*isens3
          nc   =ncl  +n0
          ncm1 =ncm1l+n0
          distn=(xcc(ncl)-xcc(ncm1l))*nxn(mbn)+ &
                (ycc(ncl)-ycc(ncm1l))*nyn(mbn)+ &
                (zcc(ncl)-zcc(ncm1l))*nzn(mbn)
          distn=dist(ncm1)+abs(distn)
          if(distn.lt.dist(nc)) then
            dist(nc) =distn
            mnpar(nc)=mbn
          else !passage au point suivant de la frontiere
            exit
          endif
          ncm1l=ncl
        enddo
!
        if(ierr.ne.0) then
          write(imp,'(/,''!!!atdist_1: probleme orientation pour '',''calcul distance'')')
          stop
        end if
      enddo
!
      return
      end
