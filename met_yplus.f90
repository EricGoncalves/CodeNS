module mod_met_yplus
implicit none
contains
      subroutine met_yplus( &
                 ncbd,ncin,s,mu,dist, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 x,y,z,nxn,nyn,nzn)
!
!***********************************************************************
!
!     ACT
!_A    Calcul de y+ au centre de la premiere maille.
!_A
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    s          : arg real(ip11,ip60 ) ; variables de calcul
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    dist       : arg real(ip12      ) ; distance a la paroi
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_L    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_L                                        contraintes visqueuses
!_L    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_L                                        contraintes visqueuses
!_L    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_L                                        contraintes visqueuses
!_L    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_L                                        contraintes visqueuses
!_L    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_L                                        contraintes visqueuses
!_L    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_L                                        contraintes visqueuses
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

!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use sortiefichier
implicit none
integer :: ncbd
integer :: ncin
double precision :: s
double precision :: dist
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: x
double precision :: y
double precision :: z
double precision :: distp
integer :: i1
integer :: i2
integer :: idf1
integer :: idf2
integer :: idfac
integer :: idm
integer :: imaxf
integer :: iminf
integer :: j1
integer :: j2
integer :: jmaxf
integer :: jminf
integer :: k1
integer :: k2
integer :: kimpl
integer :: kmaxf
integer :: kminf
integer :: l
integer :: m0b
integer :: m0n
integer :: m1
integer :: m1max
integer :: m1maxm1
integer :: m1min
integer :: m2
integer :: m2max
integer :: m2maxm1
integer :: m2min
integer :: mf
integer :: mfac
integer :: mfacn
integer :: mfl
integer :: n0c
integer :: n0n
integer :: nc
integer :: nci
integer :: ncj
integer :: nck
integer :: nfac1
integer :: nfac2
integer :: nfac3
integer :: nfac4
integer :: nfacf
integer :: nid
integer :: nijd
integer :: njd
double precision :: sxyb
double precision :: taunorm
double precision :: taupe
double precision :: utau
double precision :: utx
double precision :: utxt
double precision :: uty
double precision :: utyt
double precision :: utz
double precision :: utzt
double precision :: xcfac
double precision :: ycfac
double precision :: yplusp
double precision :: zcfac
!
!-----------------------------------------------------------------------
!
      integer dnc
      double precision nxn,nyn,nzn
      double precision mu
!
      dimension s(ip11,ip60)
      dimension x(ip12),y(ip12),z(ip12),mu(ip12),dist(ip12)
      dimension toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12)
      dimension ncbd(ip41),ncin(ip41)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
!
      kimpl=1
!
!     -------------------------------------------------------
!     SORTIES RELATIVES A DES VALEURS SUR LES PAROIS
!
      do mf=1,nbdko
!       boucle sur les parois
!
        mfl=lbdko(mf)
        l=ndlb(mfl)
!
        i1=ii1(l)
        i2=ii2(l)
        j1=jj1(l)
        j2=jj2(l)
        k1=kk1(l)
        k2=kk2(l)
        n0n=npn(l)
        n0c=npc(l)
!
        nid = id2(l)-id1(l)+1
        njd = jd2(l)-jd1(l)+1
        nijd = nid*njd
!
        nci=1
        ncj=nid
        nck=nijd
!
        m0b=mpb(mfl)
        m0n=mpn(mfl)
        iminf=iminb(mfl)
        imaxf=imaxb(mfl)
        jminf=jminb(mfl)
        jmaxf=jmaxb(mfl)
        kminf=kminb(mfl)
        kmaxf=kmaxb(mfl)
!
        if(kcaldis.eq.0) then
          write(sor2,'("!!!!met_yplus: kcaldis=0. Il faut 1 ou 2 ==>return")')
          return
        endif
!
        if(kimpl.eq.1) then
          write(sor2,987) l,mfl,iminf,imaxf,jminf,jmaxf,kminf,kmaxf
  987     format(/,'===>met_yplus: y+ pres des frontieres :', &
          /15x,24('-') &
          /1x,'zone ',i3,'  - frontiere ',i3/1x,26('-'), &
          /5x,'imin = ',i3,5x,'imax = ',i3, &
          /5x,'jmin = ',i3,5x,'jmax = ',i3, &
          /5x,'kmin = ',i3,5x,'kmax = ',i3/)
        endif
!
        m1min=1
        m2min=1
        if (iminf.eq.imaxf) then
          m1max=jmaxf-jminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          idf1=ncj
          idf2=nck
!
          if(iminf.eq.i1) then
            idfac=nci
            dnc=nci
          endif
          if(iminf.eq.i2) then
            idfac=0
            dnc=-nci
          endif
!
        elseif(jminf.eq.jmaxf) then
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          idf1=nci
          idf2=nck
!
          if(jminf.eq.j1) then
            idfac=ncj
            dnc=ncj
          endif
          if(jminf.eq.j2) then
            idfac=0
            dnc=-ncj
          endif
!
        elseif(kminf.eq.kmaxf) then
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=jmaxf-jminf+1
          m2maxm1=m2max-1
          idf1=nci
          idf2=ncj
!
          if(kminf.eq.k1) then
            idfac=nck
            dnc=nck
          endif
          if(kminf.eq.k2)then
             idfac=0
            dnc=-nck
          endif
!
        endif
!
        idm=m1max-m1min
        do m2=m2min,m2maxm1
!         boucle sur les bandes
!
!         surface bande
          sxyb=0.
!
!                 mfac   : pointeur facette tableaux toutes frontieres
!                 mfacn  : pointeur facette frontiere normale stockee
!                 nfacf  : pointeur cellule fictive
!                 nc     : pointeur cellule interne adjacente
!
!           hauteur demi maille
            write(sor2,'(" y+ et dist avec hauteur demi maille")')
            write(sor2,'("  m1  m2    mfac",t20,"x",t31,"y",t42,"z",' &
            //'t53,"utau/ai",t64,"y+_1",t75,"dist")')
!
          do m1=m1min,m1maxm1
!           boucle sur les cellules de la bande
!
            mfac=m0b+m1+(m2-1)*idm
            mfacn=m0n+m1+(m2-1)*idm
            nfac1=ncbd(mfac)-n0c+n0n+idfac
            nfac2=ncbd(mfac)-n0c+n0n+idfac+idf1
            nfac3=ncbd(mfac)-n0c+n0n+idfac+idf1+idf2
            nfac4=ncbd(mfac)-n0c+n0n+idfac+idf2
            nfacf=ncbd(mfac)
            nc=ncin(mfac)
!
            xcfac=0.25*(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))
            ycfac=0.25*(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))
            zcfac=0.25*(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))
!
!           frottement
            utx=toxx(nfacf)*nxn(mfacn)+toxy(nfacf)*nyn(mfacn)+ &
                toxz(nfacf)*nzn(mfacn)
            uty=toxy(nfacf)*nxn(mfacn)+toyy(nfacf)*nyn(mfacn)+ &
                toyz(nfacf)*nzn(mfacn)
            utz=toxz(nfacf)*nxn(mfacn)+toyz(nfacf)*nyn(mfacn)+ &
                tozz(nfacf)*nzn(mfacn)
!           projection du frottement sur la surface
            taunorm=utx*nxn(mfacn)+uty*nyn(mfacn)+utz*nzn(mfacn)
            utxt   =utx-taunorm*nxn(mfacn)
            utyt   =uty-taunorm*nyn(mfacn)
            utzt   =utz-taunorm*nzn(mfacn)
!
            taupe  =sqrt(utxt**2+utyt**2+utzt**2)
            utau   =sqrt(taupe/s(nfacf,1))
            distp  =dist(nc)
            yplusp =distp*utau*s(nc,1)/mu(nc)
!
!            hauteur demi-maille
              write(sor2,'(2i4,i8,6(1pe11.3))') &
              m1,m2,mfac,xcfac,ycfac,zcfac,utau,yplusp,distp
!
!          fin de boucle sur les cellules de la bande
           enddo
!        fin de boucle sur les bandes
         enddo
!      fin de boucle sur les parois
       enddo
!
      return
      end subroutine
end module
