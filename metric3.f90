      subroutine metric3( &
                 l,x,y,z, &
                 cvi,cvj,cvk, &
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
!
!-----------------------------------------------------------------------
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension cvi(ip21),cvj(ip21),cvk(ip21), &
                cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      indn(i,j,k)=npn(l)+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nid*njd
      inc(id,jd,kd)=id+jd*nid+kd*nid*njd
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
      n0c=npc(l)
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
!      calcul des coefficients de ponderation des schemas en maillage irregulier
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!******************************************************************************
      if(equat(3:5).eq.'2dk') then
!******************************************************************************
       do j=j1p1,j2m1
        do i=i1p1,i2m1
         n=indn(i,j,k1)
         m=n-n0c
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
!c------vecteur normal au plan interface--------------------------------
         xn=(y(n)-y(n+ncjk))*(z(n+nck)-z(n+ncj))- &
            (z(n)-z(n+ncjk))*(y(n+nck)-y(n+ncj))
         yn=(z(n)-z(n+ncjk))*(x(n+nck)-x(n+ncj))- &
            (x(n)-x(n+ncjk))*(z(n+nck)-z(n+ncj))
         zn=(x(n)-x(n+ncjk))*(y(n+nck)-y(n+ncj))- &
            (y(n)-y(n+ncjk))*(x(n+nck)-x(n+ncj))
!c-------coordonnees du point E-------------------------------------------
         tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
            (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
         xe=xp+tk*(xr-xp)
         ye=yp+tk*(yr-yp)
         ze=zp+tk*(zr-zp)
!-----------------------------------------------------------------------
!         xa1=0.5*(x(n)+x(n+nck))
!         ya1=0.5*(y(n)+y(n+nck))
!         za1=0.5*(z(n)+z(n+nck))
!         xa2=0.5*(x(n+ncj)+x(n+ncjk))
!         ya2=0.5*(y(n+ncj)+y(n+ncjk))
!         za2=0.5*(z(n+ncj)+z(n+ncjk))
!--------coefficients (attention! facteur 2 pour calcul des flux)--------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpi=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmi=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmui1(m)=2.*dpi/(dpi+dmi)
         cmui2(m)=2.*dmi/(dpi+dmi)
         cvi(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvi(m+nci)=cvi(m)*cmui2(m)
!------------------------------------------------------------------------
!      pour le flux calcule dans la direction j
!------------------------------------------------------------------------
         xp=0.125*(x(n-nck)+x(n+nci-ncj)+x(n)+x(n+ncij-ncj)+ &
                x(n+nck-ncj)+x(n+ncik-ncj)+x(n+ncjk-ncj)+x(n+ncijk-ncj))
         yp=0.125*(y(n-ncj)+y(n+nci-ncj)+y(n)+y(n+ncij-ncj)+ &
                y(n+nck-ncj)+y(n+ncik-ncj)+y(n+ncjk-ncj)+y(n+ncijk-ncj))
         zp=0.125*(z(n-ncj)+z(n+nci-ncj)+z(n)+z(n+ncij-ncj)+ &
                z(n+nck-ncj)+z(n+ncik-ncj)+z(n+ncjk-ncj)+z(n+ncijk-ncj))
!-------------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!         xa1=0.5*(x(n)+x(n+nck))
!         ya1=0.5*(y(n)+y(n+nck))
!         za1=0.5*(z(n)+z(n+nck))
!         xa2=0.5*(x(n+nci)+x(n+ncik))
!         ya2=0.5*(y(n+nci)+y(n+ncik))
!         za2=0.5*(z(n+nci)+z(n+ncik))
!-----------------------------------------------------------------------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpj=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmj=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmuj1(m)=2.*dpj/(dpj+dmj)
         cmuj2(m)=2.*dmj/(dpj+dmj)
         cvj(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvj(m+ncj)=cvj(m)*cmuj2(m)
        enddo
       enddo
!------------------------------------------------------------------------
!c----remplissage des bords
!------------------------------------------------------------------------
       do j=j1p1,j2
        m=indn(i1,j,k1)-n0c
        cmuj1(m)=cmuj1(m+nci)
        cmuj2(m)=cmuj2(m+nci)
        cvi(m)=cvi(m+nci)*cmui1(m+nci)
        cvj(m)=cvj(m+nci)
       enddo
!
       do i=i1p1,i2
        m=indn(i,j1,k1)-n0c
        cmui1(m)=cmui1(m+ncj)
        cmui2(m)=cmui2(m+ncj)
        cvj(m)=cvj(m+ncj)*cmuj1(m+ncj)
        cvi(m)=cvi(m+ncj)
       enddo
!
       mc=indn(i1,j1,k1)-n0c
       cvi(mc)=cvi(mc+ncj)
       cvj(mc)=cvj(mc+nci)
!
!******************************************************************************
      elseif(equat(3:5).eq.'2dj') then
!******************************************************************************
!
       do k=k1p1,k2m1
        do i=i1p1,i2m1
         n=indn(i,j1,k)
         m=n-n0c
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
!c------vecteur normal au plan interface--------------------------------
         xn=(y(n)-y(n+ncjk))*(z(n+nck)-z(n+ncj))- &
            (z(n)-z(n+ncjk))*(y(n+nck)-y(n+ncj))
         yn=(z(n)-z(n+ncjk))*(x(n+nck)-x(n+ncj))- &
            (x(n)-x(n+ncjk))*(z(n+nck)-z(n+ncj))
         zn=(x(n)-x(n+ncjk))*(y(n+nck)-y(n+ncj))- &
            (y(n)-y(n+ncjk))*(x(n+nck)-x(n+ncj))
!c-------coordonnees du point E-------------------------------------------
         tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
            (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
         xe=xp+tk*(xr-xp)
         ye=yp+tk*(yr-yp)
         ze=zp+tk*(zr-zp)
!-----------------------------------------------------------------------
!         xa1=0.5*(x(n)+x(n+nck))
!         ya1=0.5*(y(n)+y(n+nck))
!         za1=0.5*(z(n)+z(n+nck))
!         xa2=0.5*(x(n+ncj)+x(n+ncjk))
!         ya2=0.5*(y(n+ncj)+y(n+ncjk))
!         za2=0.5*(z(n+ncj)+z(n+ncjk))
!--------coefficients (attention! facteur 2 pour calcul des flux)--------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpi=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmi=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmui1(m)=2.*dpi/(dpi+dmi)
         cmui2(m)=2.*dmi/(dpi+dmi)
         cvi(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvi(m+nci)=cvi(m)*cmui2(m)
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
!         xa1=0.5*(x(n)+x(n+nci))
!         ya1=0.5*(y(n)+y(n+nci))
!         za1=0.5*(z(n)+z(n+nci))
!         xa2=0.5*(x(n+ncj)+x(n+ncij))
!         ya2=0.5*(y(n+ncj)+y(n+ncij))
!         za2=0.5*(z(n+ncj)+z(n+ncij))
!------------------------------------------------------------------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpk=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmk=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmuk1(m)=2.*dpk/(dpk+dmk)
         cmuk2(m)=2.*dmk/(dpk+dmk)
         cvk(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvk(m+nck)=cvk(m)*cmuk2(m)
        enddo
       enddo
!------------------------------------------------------------------------
!c----remplissage des bords
!------------------------------------------------------------------------
       do k=k1p1,k2
        m=indn(i1,j1,k)-n0c
        cmuk1(m)=cmuk1(m+nci)
        cmuk2(m)=cmuk2(m+nci)
        cvi(m)=cvi(m+nci)*cmui1(m+nci)
        cvk(m)=cvk(m+nci)
       enddo
!
       do i=i1p1,i2
        m=indn(i,j1,k1)-n0c
        cmui1(m)=cmui1(m+nck)
        cmui2(m)=cmui2(m+nck)
        cvk(m)=cvk(m+nck)*cmuk1(m+nck)
        cvi(m)=cvi(m+nck)
       enddo
!
       mc=indn(i1,j1,k1)-n0c
       cvi(mc)=cvi(mc+nck)
       cvk(mc)=cvk(mc+nci)
!
!******************************************************************************
      elseif(equat(3:5).eq.'2di') then
!******************************************************************************
!
       do k=k1p1,k2m1
        do j=j1p1,j2m1
         n=indn(i1,j,k)
         m=n-n0c
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
!-------------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!         xa1=0.5*(x(n)+x(n+nck))
!         ya1=0.5*(y(n)+y(n+nck))
!         za1=0.5*(z(n)+z(n+nck))
!         xa2=0.5*(x(n+nci)+x(n+ncik))
!         ya2=0.5*(y(n+nci)+y(n+ncik))
!         za2=0.5*(z(n+nci)+z(n+ncik))
!-----------------------------------------------------------------------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpj=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmj=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmuj1(m)=2.*dpj/(dpj+dmj)
         cmuj2(m)=2.*dmj/(dpj+dmj)
         cvj(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvj(m+ncj)=cvj(m)*cmuj2(m)
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
!         xa1=0.5*(x(n)+x(n+nci))
!         ya1=0.5*(y(n)+y(n+nci))
!         za1=0.5*(z(n)+z(n+nci))
!         xa2=0.5*(x(n+ncj)+x(n+ncij))
!         ya2=0.5*(y(n+ncj)+y(n+ncij))
!         za2=0.5*(z(n+ncj)+z(n+ncij))
!------------------------------------------------------------------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpk=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmk=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmuk1(m)=2.*dpk/(dpk+dmk)
         cmuk2(m)=2.*dmk/(dpk+dmk)
         cvk(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvk(m+nck)=cvk(m)*cmuk2(m)
        enddo
       enddo
!------------------------------------------------------------------------
!c----remplissage des bords
!------------------------------------------------------------------------
       do j=j1p1,j2
        m=indn(i1,j,k1)-n0c
        cmuj1(m)=cmuj1(m+nck)
        cmuj2(m)=cmuj2(m+nck)
        cvk(m)=cvk(m+nck)*cmuk1(m+nck)
        cvj(m)=cvj(m+nck)
       enddo
!
       do k=k1p1,k2
        m=indn(i1,j1,k)-n0c
        cmuk1(m)=cmuk1(m+ncj)
        cmuk2(m)=cmuk2(m+ncj)
        cvj(m)=cvj(m+ncj)*cmuj1(m+ncj)
        cvk(m)=cvk(m+ncj)
       enddo
!
       mc=indn(i1,j1,k1)-n0c
       cvj(mc)=cvj(mc+nck)
       cvk(mc)=cvk(mc+ncj)
!
!******************************************************************************
      elseif(equat(3:4).eq.'3d') then
!******************************************************************************
!
      do k=k1p1,k2m1
       do j=j1p1,j2m1
        do i=i1p1,i2m1
         n=indn(i,j,k)
         m=n-n0c
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
!c------vecteur normal au plan interface--------------------------------
         xn=(y(n)-y(n+ncjk))*(z(n+nck)-z(n+ncj))- &
            (z(n)-z(n+ncjk))*(y(n+nck)-y(n+ncj))
         yn=(z(n)-z(n+ncjk))*(x(n+nck)-x(n+ncj))- &
            (x(n)-x(n+ncjk))*(z(n+nck)-z(n+ncj))
         zn=(x(n)-x(n+ncjk))*(y(n+nck)-y(n+ncj))- &
            (y(n)-y(n+ncjk))*(x(n+nck)-x(n+ncj))
!c-------coordonnees du point E-------------------------------------------
         tk=(xn*(x(n)-xp)+yn*(y(n)-yp)+zn*(z(n)-zp))/ &
            (xn*(xr-xp)+yn*(yr-yp)+zn*(zr-zp))
         xe=xp+tk*(xr-xp)
         ye=yp+tk*(yr-yp)
         ze=zp+tk*(zr-zp)
!-----------------------------------------------------------------------
!         xa1=0.5*(x(n)+x(n+nck))
!         ya1=0.5*(y(n)+y(n+nck))
!         za1=0.5*(z(n)+z(n+nck))
!         xa2=0.5*(x(n+ncj)+x(n+ncjk))
!         ya2=0.5*(y(n+ncj)+y(n+ncjk))
!         za2=0.5*(z(n+ncj)+z(n+ncjk))
!--------coefficients (attention! facteur 2 pour calcul des flux)--------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpi=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmi=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmui1(m)=2.*dpi/(dpi+dmi)
         cmui2(m)=2.*dmi/(dpi+dmi)
         cvi(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvi(m+nci)=cvi(m)*cmui2(m)
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
!----------------------------------------------------------------------------
!         xa1=0.5*(x(n)+x(n+nck))
!         ya1=0.5*(y(n)+y(n+nck))
!         za1=0.5*(z(n)+z(n+nck))
!         xa2=0.5*(x(n+nci)+x(n+ncik))
!         ya2=0.5*(y(n+nci)+y(n+ncik))
!         za2=0.5*(z(n+nci)+z(n+ncik))
!----------------------------------------------------------------------------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpj=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmj=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmuj1(m)=2.*dpj/(dpj+dmj)
         cmuj2(m)=2.*dmj/(dpj+dmj)
         cvj(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvj(m+ncj)=cvj(m)*cmuj2(m)
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
!         xa1=0.5*(x(n)+x(n+nci))
!         ya1=0.5*(y(n)+y(n+nci))
!         za1=0.5*(z(n)+z(n+nci))
!         xa2=0.5*(x(n+ncj)+x(n+ncij))
!         ya2=0.5*(y(n+ncj)+y(n+ncij))
!         za2=0.5*(z(n+ncj)+z(n+ncij))
!------------------------------------------------------------------
!         d2=sqrt((xa1-xe)**2+(ya1-ye)**2+(za1-ze)**2)
!         d1=sqrt((xa2-xe)**2+(ya2-ye)**2+(za2-ze)**2)
!         d3=d1+d2
!         d1=d1/d3
!         d2=d2/d3
         dpk=sqrt((xp-xe)**2+(yp-ye)**2+(zp-ze)**2)
         dmk=sqrt((xr-xe)**2+(yr-ye)**2+(zr-ze)**2)
         cmuk1(m)=2.*dpk/(dpk+dmk)
         cmuk2(m)=2.*dmk/(dpk+dmk)
         cvk(m)=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
         cvk(m+nck)=cvk(m)*cmuk2(m)
        enddo
       enddo
      enddo
!------------------------------------------------------------------------
!c----remplissage des bords
!------------------------------------------------------------------------
      do j=j1p1,j2
       do k=k1p1,k2
        m=indn(i1,j,k)-n0c
        cmuj1(m)=cmuj1(m+nci)
        cmuj2(m)=cmuj2(m+nci)
        cmuk1(m)=cmuk1(m+nci)
        cmuk2(m)=cmuk2(m+nci)
        cvi(m)=cvi(m+nci)*cmui1(m+nci)
        cvj(m)=cvj(m+nci)
        cvk(m)=cvk(m+nci)
       enddo
      enddo
!
      do i=i1p1,i2
       do j=j1p1,j2
        m=indn(i,j,k1)-n0c
        cmui1(m)=cmui1(m+nck)
        cmui2(m)=cmui2(m+nck)
        cmuj1(m)=cmuj1(m+nck)
        cmuj2(m)=cmuj2(m+nck)
        cvk(m)=cvk(m+nck)*cmuk1(m+nck)
        cvi(m)=cvi(m+nck)
        cvj(m)=cvi(m+nck)
       enddo
      enddo
!
      do i=i1p1,i2
       do k=k1p1,k2
        m=indn(i,j1,k)-n0c
        cmui1(m)=cmui1(m+ncj)
        cmui2(m)=cmui2(m+ncj)
        cmuk1(m)=cmuk1(m+ncj)
        cmuk2(m)=cmuk2(m+ncj)
        cvj(m)=cvj(m+ncj)*cmuj1(m+ncj)
        cvi(m)=cvi(m+ncj)
        cvk(m)=cvk(m+ncj)
       enddo
      enddo
!
      do i=i1p1,i2
       m=indn(i,j1,k1)-n0c
       cmui1(m)=cmui1(m+ncj)
       cmui2(m)=cmui2(m+ncj)
       cvj(m)=cvj(m+ncj)*cmuj1(m+ncj)
       cvk(m)=cvk(m+nck)*cmuk1(m+nck)
       cvi(m)=cvi(m+ncj)
      enddo
!
      do j=j1p1,j2
       m=indn(i1,j,k1)-n0c
       cmuj1(m)=cmuj1(m+nci)
       cmuj2(m)=cmuj2(m+nci)
       cvi(m)=cvi(m+nci)*cmui1(m+nci)
       cvk(m)=cvk(m+nck)*cmuk1(m+nck)
       cvj(m)=cvj(m+nci)
      enddo
!
      do k=k1p1,k2
       m=indn(i1,j1,k)-n0c
       cmuk1(m)=cmuk1(m+ncj)
       cmuk2(m)=cmuk2(m+ncj)
       cvi(m)=cvi(m+nci)*cmui1(m+nci)
       cvj(m)=cvj(m+ncj)*cmuj1(m+ncj)
       cvk(m)=cvk(m+ncj)
      enddo
!
       mc=indn(i1,j1,k1)-n0c
       cvi(mc)=cvi(mc+ncj)
       cvj(mc)=cvj(mc+nci)
       cvk(mc)=cvk(mc+ncj)
!******************************************************************
      endif
!
      return
      end
