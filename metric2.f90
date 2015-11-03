module mod_metric2
  implicit none
contains
  subroutine metric2( &
       l,x,y,z, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C : decembre 2001 - Eric GONCALVES / SINUMEF
!
!     ACT
!_A     Coefficients de ponderation pour maillage irregulier.
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
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
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    implicit none
    integer          ::     i,   i1, i1p1,   i2, i2m1
    integer          ::    id,    j,   j1, j1p1,   j2
    integer          ::  j2m1,   jd,    k,   k1, k1p1
    integer          ::    k2, k2m1,   kd,    l    
    integer          ::     n,  nci, ncij,ncijk
    integer          ::  ncik,  ncj, ncjk,  nck,  nid
    integer          ::   njd
    double precision :: cmui1,cmui2,cmuj1,cmuj2,cmuk1
    double precision :: cmuk2,  dmi,  dmj,  dmk,  dpi
    double precision ::   dpj,  dpk,   tk,    x,   xe
    double precision ::    xn,   xp,   xr,    y,   ye
    double precision ::    yn,   yp,   yr,    z,   ze
    double precision ::    zn,   zp,   zr
!
!-----------------------------------------------------------------------
!
    dimension x(ip21),y(ip21),z(ip21)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
              cmuk1(ip21),cmuk2(ip21)
!
    nid=id2(l)-id1(l)+1
    njd=jd2(l)-jd1(l)+1
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
    ncij= inc(1,1,0)
    ncik= inc(1,0,1)
    ncjk= inc(0,1,1)
    ncijk= inc(1,1,1)
!
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
    i1p1=i1+1
    j1p1=j1+1
    k1p1=k1+1
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      calcul des coefficients bi1 (bi2) multiplicateurs a la derivee 1ere (2eme)
!      dans les directions transversales au flux calcule
!      calcul des coefficients  dp,dm,dc utile a l'approximation de la
!      derivee 2eme dans la direction du flux calcule
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!******************************************************************************
    if(equat(3:5).eq.'2dk') then
!******************************************************************************
       do j=j1p1,j2m1
          do i=i1p1,i2m1
             n=indn(i,j,k1)
!-----------------------------------------------------------------------
!             pour le flux calcule dans la direction i
!-----------------------------------------------------------------------
!       points P,Q,R centres cellules
!       points E,F points interfaces cellules
!
!    (n-nci)     (n)     (n+nci)
!      P --- E ---R ------ Q
!
!      distance EP -> dp
!      distance ER -> dm
!      distance RP -> cv
!
             xp=0.125*(x(n-nci)+x(n)+x(n+ncj-nci)+x(n+ncij-nci)+ &
                       x(n+nck-nci)+x(n+ncik-nci)+x(n+ncjk-nci)+x(n+ncijk-nci))
             yp=0.125*(y(n-nci)+y(n)+y(n+ncj-nci)+y(n+ncij-nci)+ &
                       y(n+nck-nci)+y(n+ncik-nci)+y(n+ncjk-nci)+y(n+ncijk-nci))
             zp=0.125*(z(n-nci)+z(n)+z(n+ncj-nci)+z(n+ncij-nci)+ &
                       z(n+nck-nci)+z(n+ncik-nci)+z(n+ncjk-nci)+z(n+ncijk-nci))
             xr=0.125*(x(n)+x(n+nci)+x(n+ncj)+x(n+ncij)+ &
                       x(n+nck)+x(n+ncik)+x(n+ncjk)+x(n+ncijk))
             yr=0.125*(y(n)+y(n+nci)+y(n+ncj)+y(n+ncij)+ &
                       y(n+nck)+y(n+ncik)+y(n+ncjk)+y(n+ncijk))
             zr=0.125*(z(n)+z(n+nci)+z(n+ncj)+z(n+ncij)+ &
                       z(n+nck)+z(n+ncik)+z(n+ncjk)+z(n+ncijk))
!-------vecteur normal au plan interface--------------------------------
             xn=(y(n)-y(n+ncjk))*(z(n+nck)-z(n+ncj))- &
                (z(n)-z(n+ncjk))*(y(n+nck)-y(n+ncj))
             yn=(z(n)-z(n+ncjk))*(x(n+nck)-x(n+ncj))- &
                (x(n)-x(n+ncjk))*(z(n+nck)-z(n+ncj))
             zn=(x(n)-x(n+ncjk))*(y(n+nck)-y(n+ncj))- &
                (y(n)-y(n+ncjk))*(x(n+nck)-x(n+ncj))
!--------coordonnees du point E-------------------------------------------
             tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                  (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
             xe=xp+tk*(xr-xp)
             ye=yp+tk*(yr-yp)
             ze=zp+tk*(zr-zp)
             dpi=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
             dmi=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
!--------pour operateur mu pondere------------------------------
!        attention! facteur 2 pour calcul des flux
             cmui1(n)=2.*dpi/(dpi+dmi)
             cmui2(n)=2.*dmi/(dpi+dmi)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction j
!------------------------------------------------------------------------
             xp=0.125*(x(n-ncj)+x(n+nci-ncj)+x(n)+x(n+ncij-ncj)+ &
                       x(n+nck-ncj)+x(n+ncik-ncj)+x(n+ncjk-ncj)+x(n+ncijk-ncj))
             yp=0.125*(y(n-ncj)+y(n+nci-ncj)+y(n)+y(n+ncij-ncj)+ &
                       y(n+nck-ncj)+y(n+ncik-ncj)+y(n+ncjk-ncj)+y(n+ncijk-ncj))
             zp=0.125*(z(n-ncj)+z(n+nci-ncj)+z(n)+z(n+ncij-ncj)+ &
                       z(n+nck-ncj)+z(n+ncik-ncj)+z(n+ncjk-ncj)+z(n+ncijk-ncj))
!--------------------------------------------------------------------
             xn=(y(n)-y(n+ncik))*(z(n+nck)-z(n+nci))- &
                (z(n)-z(n+ncik))*(y(n+nck)-y(n+nci))
             yn=(z(n)-z(n+ncik))*(x(n+nck)-x(n+nci))- &
                (x(n)-x(n+ncik))*(z(n+nck)-z(n+nci))
             zn=(x(n)-x(n+ncik))*(y(n+nck)-y(n+nci))- &
                (y(n)-y(n+ncik))*(x(n+nck)-x(n+nci))
             tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
             xe=xp+tk*(xr-xp)
             ye=yp+tk*(yr-yp)
             ze=zp+tk*(zr-zp)
             dpj=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
             dmj=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
             cmuj1(n)=2.*dpj/(dpj+dmj)
             cmuj2(n)=2.*dmj/(dpj+dmj)
          enddo
       enddo
!------------------------------------------------------------------------
!-----remplissage des bords
!------------------------------------------------------------------------
       do j=j1p1,j2m1
          n=indn(i1,j,k1)
          cmuj1(n)=cmuj1(n+nci)
          cmuj2(n)=cmuj2(n+nci)
       enddo
!
       do i=i1p1,i2m1
          n=indn(i,j1,k1)
          cmui1(n)=cmui1(n+ncj)
          cmui2(n)=cmui2(n+ncj)
       enddo
!
!******************************************************************************
    elseif(equat(3:5).eq.'2dj') then
!******************************************************************************
!
       do k=k1p1,k2m1
          do i=i1p1,i2m1
             n=indn(i,j1,k)
!-----------------------------------------------------------------------
!             pour le flux calcule dans la direction i
!-----------------------------------------------------------------------
             xp=0.125*(x(n-nci)+x(n)+x(n+ncj-nci)+x(n+ncij-nci)+ &
                       x(n+nck-nci)+x(n+ncik-nci)+x(n+ncjk-nci)+x(n+ncijk-nci))
             yp=0.125*(y(n-nci)+y(n)+y(n+ncj-nci)+y(n+ncij-nci)+ &
                       y(n+nck-nci)+y(n+ncik-nci)+y(n+ncjk-nci)+y(n+ncijk-nci))
             zp=0.125*(z(n-nci)+z(n)+z(n+ncj-nci)+z(n+ncij-nci)+ &
                       z(n+nck-nci)+z(n+ncik-nci)+z(n+ncjk-nci)+z(n+ncijk-nci))
             xr=0.125*(x(n)+x(n+nci)+x(n+ncj)+x(n+ncij)+ &
                       x(n+nck)+x(n+ncik)+x(n+ncjk)+x(n+ncijk))
             yr=0.125*(y(n)+y(n+nci)+y(n+ncj)+y(n+ncij)+ &
                       y(n+nck)+y(n+ncik)+y(n+ncjk)+y(n+ncijk))
             zr=0.125*(z(n)+z(n+nci)+z(n+ncj)+z(n+ncij)+ &
                       z(n+nck)+z(n+ncik)+z(n+ncjk)+z(n+ncijk))
!-------vecteur normal au plan interface--------------------------------
             xn=(y(n)-y(n+ncjk))*(z(n+nck)-z(n+ncj))- &
                (z(n)-z(n+ncjk))*(y(n+nck)-y(n+ncj))
             yn=(z(n)-z(n+ncjk))*(x(n+nck)-x(n+ncj))- &
                (x(n)-x(n+ncjk))*(z(n+nck)-z(n+ncj))
             zn=(x(n)-x(n+ncjk))*(y(n+nck)-y(n+ncj))- &
                (y(n)-y(n+ncjk))*(x(n+nck)-x(n+ncj))
!--------coordonnees du point E-------------------------------------------
             tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
             xe=xp+tk*(xr-xp)
             ye=yp+tk*(yr-yp)
             ze=zp+tk*(zr-zp)
             dpi=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
             dmi=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
!--------pour operateur mu pondere------------------------------
!        attention! facteur 2 pour calcul des flux
             cmui1(n)=2.*dpi/(dpi+dmi)
             cmui2(n)=2.*dmi/(dpi+dmi)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction k
!------------------------------------------------------------------------
             xp=0.125*(x(n-nck)+x(n+nci-nck)+x(n)+x(n+nci)+ &
                       x(n-nck+ncj)+x(n+ncj)+x(n+ncij-nck)+x(n+ncij))
             yp=0.125*(y(n-nck)+y(n+nci-nck)+y(n)+y(n+nci)+ &
                       y(n-nck+ncj)+y(n+ncj)+y(n+ncij-nck)+y(n+ncij))
             zp=0.125*(z(n-nck)+z(n+nci-nck)+z(n)+z(n+nci)+ &
                       z(n-nck+ncj)+z(n+ncj)+z(n+ncij-nck)+z(n+ncij))
!-----------------------------------------------------------------------
             xn=(y(n)-y(n+ncij))*(z(n+ncj)-z(n+nci))- &
                (z(n)-z(n+ncij))*(y(n+ncj)-y(n+nci))
             yn=(z(n)-z(n+ncij))*(x(n+ncj)-x(n+nci))- &
                (x(n)-x(n+ncij))*(z(n+ncj)-z(n+nci))
             zn=(x(n)-x(n+ncij))*(y(n+ncj)-y(n+nci))- &
                (y(n)-y(n+ncij))*(x(n+ncj)-x(n+nci))
             tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
             xe=xp+tk*(xr-xp)
             ye=yp+tk*(yr-yp)
             ze=zp+tk*(zr-zp)
!------------------------------------------------------------------
             dpk=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
             dmk=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
             cmuk1(n)=2.*dpk/(dpk+dmk)
             cmuk2(n)=2.*dmk/(dpk+dmk)
          enddo
       enddo
!------------------------------------------------------------------------
!-----remplissage des bords
!------------------------------------------------------------------------
       do k=k1p1,k2m1
          n=indn(i1,j1,k)
          cmuk1(n)=cmuk1(n+nci)
          cmuk2(n)=cmuk2(n+nci)
       enddo
!
       do i=i1p1,i2m1
          n=indn(i,j1,k1)
          cmui1(n)=cmui1(n+nck)
          cmui2(n)=cmui2(n+nck)
       enddo
!
!******************************************************************************
    elseif(equat(3:5).eq.'2di') then
!******************************************************************************
!
       do k=k1p1,k2m1
          do j=j1p1,j2m1
             n=indn(i1,j,k)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction j
!------------------------------------------------------------------------
             xp=0.125*(x(n-ncj)+x(n+nci-ncj)+x(n)+x(n+ncij-ncj)+ &
                       x(n+nck-ncj)+x(n+ncik-ncj)+x(n+ncjk-ncj)+x(n+ncijk-ncj))
             yp=0.125*(y(n-ncj)+y(n+nci-ncj)+y(n)+y(n+ncij-ncj)+ &
                       y(n+nck-ncj)+y(n+ncik-ncj)+y(n+ncjk-ncj)+y(n+ncijk-ncj))
             zp=0.125*(z(n-ncj)+z(n+nci-ncj)+z(n)+z(n+ncij-ncj)+ &
                       z(n+nck-ncj)+z(n+ncik-ncj)+z(n+ncjk-ncj)+z(n+ncijk-ncj))
             xr=0.125*(x(n)+x(n+nci)+x(n+ncj)+x(n+ncij)+ &
                       x(n+nck)+x(n+ncik)+x(n+ncjk)+x(n+ncijk))
             yr=0.125*(y(n)+y(n+nci)+y(n+ncj)+y(n+ncij)+ &
                       y(n+nck)+y(n+ncik)+y(n+ncjk)+y(n+ncijk))
             zr=0.125*(z(n)+z(n+nci)+z(n+ncj)+z(n+ncij)+ &
                       z(n+nck)+z(n+ncik)+z(n+ncjk)+z(n+ncijk))
!--------------------------------------------------------------------
             xn=(y(n)-y(n+ncik))*(z(n+nck)-z(n+nci))- &
                (z(n)-z(n+ncik))*(y(n+nck)-y(n+nci))
             yn=(z(n)-z(n+ncik))*(x(n+nck)-x(n+nci))- &
                (x(n)-x(n+ncik))*(z(n+nck)-z(n+nci))
             zn=(x(n)-x(n+ncik))*(y(n+nck)-y(n+nci))- &
                (y(n)-y(n+ncik))*(x(n+nck)-x(n+nci))
             tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
             xe=xp+tk*(xr-xp)
             ye=yp+tk*(yr-yp)
             ze=zp+tk*(zr-zp)
             dpj=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
             dmj=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
             cmuj1(n)=2.*dpj/(dpj+dmj)
             cmuj2(n)=2.*dmj/(dpj+dmj)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction k
!------------------------------------------------------------------------
             xp=0.125*(x(n-nck)+x(n+nci-nck)+x(n)+x(n+nci)+ &
                       x(n-nck+ncj)+x(n+ncj)+x(n+ncij-nck)+x(n+ncij))
             yp=0.125*(y(n-nck)+y(n+nci-nck)+y(n)+y(n+nci)+ &
                       y(n-nck+ncj)+y(n+ncj)+y(n+ncij-nck)+y(n+ncij))
             zp=0.125*(z(n-nck)+z(n+nci-nck)+z(n)+z(n+nci)+ &
                       z(n-nck+ncj)+z(n+ncj)+z(n+ncij-nck)+z(n+ncij))
!-----------------------------------------------------------------------
             xn=(y(n)-y(n+ncij))*(z(n+ncj)-z(n+nci))- &
                (z(n)-z(n+ncij))*(y(n+ncj)-y(n+nci))
             yn=(z(n)-z(n+ncij))*(x(n+ncj)-x(n+nci))- &
                (x(n)-x(n+ncij))*(z(n+ncj)-z(n+nci))
             zn=(x(n)-x(n+ncij))*(y(n+ncj)-y(n+nci))- &
                (y(n)-y(n+ncij))*(x(n+ncj)-x(n+nci))
             tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
             xe=xp+tk*(xr-xp)
             ye=yp+tk*(yr-yp)
             ze=zp+tk*(zr-zp)
!------------------------------------------------------------------
             dpk=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
             dmk=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
             cmuk1(n)=2.*dpk/(dpk+dmk)
             cmuk2(n)=2.*dmk/(dpk+dmk)
          enddo
       enddo
!------------------------------------------------------------------------
!-----remplissage des bords
!------------------------------------------------------------------------
       do j=j1p1,j2m1
          n=indn(i1,j,k1)
          cmuj1(n)=cmuj1(n+nck)
          cmuj2(n)=cmuj2(n+nck)
       enddo
!
       do k=k1p1,k2m1
          n=indn(i1,j1,k)
          cmuk1(n)=cmuk1(n+ncj)
          cmuk2(n)=cmuk2(n+ncj)
       enddo
!
!******************************************************************************
    elseif(equat(3:4).eq.'3d') then
!******************************************************************************
!
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
               n=indn(i,j,k)
!-----------------------------------------------------------------------
!             pour le flux calcule dans la direction i
!-----------------------------------------------------------------------
!       points P,Q,R centres cellules
!       points E,F points interfaces cellules
!
!    (n-nci)     (n)     (n+nci)
!      P --- E ---R ------ Q
!
!      distance EP -> dp
!      distance ER -> dm
!      distance RP -> cv
!
                xp=0.125*(x(n-nci)+x(n)+x(n+ncj-nci)+x(n+ncij-nci)+ &
                          x(n+nck-nci)+x(n+ncik-nci)+x(n+ncjk-nci)+x(n+ncijk-nci))
                yp=0.125*(y(n-nci)+y(n)+y(n+ncj-nci)+y(n+ncij-nci)+ &
                          y(n+nck-nci)+y(n+ncik-nci)+y(n+ncjk-nci)+y(n+ncijk-nci))
                zp=0.125*(z(n-nci)+z(n)+z(n+ncj-nci)+z(n+ncij-nci)+ &
                          z(n+nck-nci)+z(n+ncik-nci)+z(n+ncjk-nci)+z(n+ncijk-nci))
                xr=0.125*(x(n)+x(n+nci)+x(n+ncj)+x(n+ncij)+ &
                          x(n+nck)+x(n+ncik)+x(n+ncjk)+x(n+ncijk))
                yr=0.125*(y(n)+y(n+nci)+y(n+ncj)+y(n+ncij)+ &
                          y(n+nck)+y(n+ncik)+y(n+ncjk)+y(n+ncijk))
                zr=0.125*(z(n)+z(n+nci)+z(n+ncj)+z(n+ncij)+ &
                          z(n+nck)+z(n+ncik)+z(n+ncjk)+z(n+ncijk))
!-------vecteur normal au plan interface--------------------------------
                xn=(y(n)-y(n+ncjk))*(z(n+nck)-z(n+ncj))- &
                   (z(n)-z(n+ncjk))*(y(n+nck)-y(n+ncj))
                yn=(z(n)-z(n+ncjk))*(x(n+nck)-x(n+ncj))- &
                   (x(n)-x(n+ncjk))*(z(n+nck)-z(n+ncj))
                zn=(x(n)-x(n+ncjk))*(y(n+nck)-y(n+ncj))- &
                   (y(n)-y(n+ncjk))*(x(n+nck)-x(n+ncj))
!--------coordonnees du point E-------------------------------------------
                tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                   (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
                xe=xp+tk*(xr-xp)
                ye=yp+tk*(yr-yp)
                ze=zp+tk*(zr-zp)
                dpi=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
                dmi=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
                cmui1(n)=2.*dpi/(dpi+dmi)
                cmui2(n)=2.*dmi/(dpi+dmi)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction j
!------------------------------------------------------------------------
                xp=0.125*(x(n-nck)+x(n+nci-ncj)+x(n)+x(n+ncij-ncj)+ &
                          x(n+nck-ncj)+x(n+ncik-ncj)+x(n+ncjk-ncj)+x(n+ncijk-ncj))
                yp=0.125*(y(n-ncj)+y(n+nci-ncj)+y(n)+y(n+ncij-ncj)+ &
                          y(n+nck-ncj)+y(n+ncik-ncj)+y(n+ncjk-ncj)+y(n+ncijk-ncj))
                zp=0.125*(z(n-ncj)+z(n+nci-ncj)+z(n)+z(n+ncij-ncj)+ &
                          z(n+nck-ncj)+z(n+ncik-ncj)+z(n+ncjk-ncj)+z(n+ncijk-ncj))
!--------------------------------------------------------------------
                xn=(y(n)-y(n+ncik))*(z(n+nck)-z(n+nci))- &
                   (z(n)-z(n+ncik))*(y(n+nck)-y(n+nci))
                yn=(z(n)-z(n+ncik))*(x(n+nck)-x(n+nci))- &
                   (x(n)-x(n+ncik))*(z(n+nck)-z(n+nci))
                zn=(x(n)-x(n+ncik))*(y(n+nck)-y(n+nci))- &
                   (y(n)-y(n+ncik))*(x(n+nck)-x(n+nci))
                tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                   (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
                xe=xp+tk*(xr-xp)
                ye=yp+tk*(yr-yp)
                ze=zp+tk*(zr-zp)
                dpj=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
                dmj=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
                cmuj1(n)=2.*dpj/(dpj+dmj)
                cmuj2(n)=2.*dmj/(dpj+dmj)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction k
!------------------------------------------------------------------------
                xp=0.125*(x(n-nck)+x(n+nci-nck)+x(n)+x(n+nci)+ &
                          x(n-nck+ncj)+x(n+ncj)+x(n+ncij-nck)+x(n+ncij))
                yp=0.125*(y(n-nck)+y(n+nci-nck)+y(n)+y(n+nci)+ &
                          y(n-nck+ncj)+y(n+ncj)+y(n+ncij-nck)+y(n+ncij))
                zp=0.125*(z(n-nck)+z(n+nci-nck)+z(n)+z(n+nci)+ &
                          z(n-nck+ncj)+z(n+ncj)+z(n+ncij-nck)+z(n+ncij))
!-----------------------------------------------------------------------
                xn=(y(n)-y(n+ncij))*(z(n+ncj)-z(n+nci))- &
                   (z(n)-z(n+ncij))*(y(n+ncj)-y(n+nci))
                yn=(z(n)-z(n+ncij))*(x(n+ncj)-x(n+nci))- &
                   (x(n)-x(n+ncij))*(z(n+ncj)-z(n+nci))
                zn=(x(n)-x(n+ncij))*(y(n+ncj)-y(n+nci))- &
                   (y(n)-y(n+ncij))*(x(n+ncj)-x(n+nci))
                tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
                   (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
                xe=xp+tk*(xr-xp)
                ye=yp+tk*(yr-yp)
                ze=zp+tk*(zr-zp)
!------------------------------------------------------------------
                dpk=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
                dmk=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
                cmuk1(n)=2.*dpk/(dpk+dmk)
                cmuk2(n)=2.*dmk/(dpk+dmk)
             enddo
          enddo
       enddo
!------------------------------------------------------------------------
!-----remplissage des bords
!------------------------------------------------------------------------
       do j=j1p1,j2m1
          do k=k1p1,k2m1
             n=indn(i1,j,k)
             cmuj1(n)=cmuj1(n+nci)
             cmuj2(n)=cmuj2(n+nci)
             cmuk1(n)=cmuk1(n+nci)
             cmuk2(n)=cmuk2(n+nci)
          enddo
       enddo
!
       do i=i1p1,i2m1
          do j=j1p1,j2m1
             n=indn(i,j,k1)
             cmui1(n)=cmui1(n+nck)
             cmui2(n)=cmui2(n+nck)
             cmuj1(n)=cmuj1(n+nck)
             cmuj2(n)=cmuj2(n+nck)
          enddo
       enddo
!
       do i=i1p1,i2m1
          do k=k1p1,k2m1
             n=indn(i,j1,k)
             cmui1(n)=cmui1(n+ncj)
             cmui2(n)=cmui2(n+ncj)
             cmuk1(n)=cmuk1(n+ncj)
             cmuk2(n)=cmuk2(n+ncj)
          enddo
       enddo
!
       do i=i1p1,i2m1
          n=indn(i,j1,k1)
          cmui1(n)=cmui1(n+ncj)
          cmui2(n)=cmui2(n+ncj)
       enddo
!
       do j=j1p1,j2m1
          n=indn(i1,j,k1)
          cmuj1(n)=cmuj1(n+nci)
          cmuj2(n)=cmuj2(n+nci)
       enddo
!
       do k=k1p1,k2m1
          n=indn(i1,j1,k)
          cmuk1(n)=cmuk1(n+ncj)
          cmuk2(n)=cmuk2(n+ncj)
       enddo
!******************************************************************
    endif
!
    return
  contains
    function    indn(i,j,k)
      implicit none
      integer          ::    i,indn,   j,   k
      indn=npn(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nid*njd
    end function indn
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nid*njd
    end function inc
  end subroutine metric2
end module mod_metric2
