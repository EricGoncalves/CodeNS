module mod_implimf_3d
  implicit none
contains
  subroutine implimf_3d( &
       lm,u,dt,v,d,ff, &
       mu,mut, &
       equat,lmx, &
       sn,lgsnlt, &
       vol,dtpas,ityprk, &
       dfxx,dfyy,dfzz,dfxy,dfxz,dfyz,dfex,dfey,dfez,rv, &
       ps,cson)
!
!***********************************************************************
!
!_DA  DATE_C : mars 2002 - Eric Goncalves / SINUMEF
!
!     ACT
!_A    Phase implicite sans matrice avec relaxation
!_A    au moyen d'une methode de Jacobi par points.
!
!***********************************************************************
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use schemanum
    implicit none
    integer          ::     i1,  i1m1,    i2,  i2m1
    integer          ::   ind1,  ind2,ityprk,     j
    integer          ::     j1,  j1m1,    j2,  j2m1
    integer          ::      k,    k1,  k1m1,    k2,  k2m1
    integer          ::   kdir,lgsnlt,    lm,   lmx
    integer          ::     ls,     m,     n,   n0c,   nci
    integer          ::    ncj,   nck,   nid,  nijd,  ninc
    integer          ::    njd
    double precision ::                   cc,                cnds,          cson(ip11),        d(ip11,ip60),          dfex(ip00)
    double precision ::           dfey(ip00),          dfez(ip00),          dfxx(ip00),          dfxy(ip00),          dfxz(ip00)
    double precision ::           dfyy(ip00),          dfyz(ip00),          dfzz(ip00),            dt(ip11),               dtpas
    double precision ::                 fact,                 fex,                 fey,                 fez,       ff(ip11,ip60)
    double precision ::                  fxx,                 fxy,                 fxz,                 fyy,                 fyz
    double precision ::                  fzz,            mu(ip12),           mut(ip12),                pres,            ps(ip11)
    double precision ::             rv(ip00),sn(lgsnlt,nind,ndir),                 tn1,                 tn2,                 tn3
    double precision ::                  tn4,                 tn5,        u(ip11,ip60),                  ui,                  uu
    double precision ::         v(ip11,ip60),                  vi,                  vn,           vol(ip11),                  vv
    double precision ::                   wi,                 wi1,                 wi2,                 wi3,                 wi4
    double precision ::                  wi5,                  ww
    double precision,allocatable :: coefdiag(:), coefe(:,:),    d2w1(:),    d2w2(:),    d2w3(:)
    double precision,allocatable ::     d2w4(:),    d2w5(:)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!
    ALLOCATE(coefe(ndir,ip00))
    ALLOCATE(coefdiag(ip00),d2w1(ip00),d2w2(ip00),d2w3(ip00),d2w4(ip00),d2w5(ip00))

    n0c=npc(lm)
    i1=ii1(lm)
    i2=ii2(lm)
    j1=jj1(lm)
    j2=jj2(lm)
    k1=kk1(lm)
    k2=kk2(lm)
!
    nid = id2(lm)-id1(lm)+1
    njd = jd2(lm)-jd1(lm)+1
    nijd = nid*njd
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)

!   constante instationnaire DTS
    fact=1.5
!   fact=11./6.  !ordre 3

!-----initalisation--------------------------------
!
    ind1 = indc(i1m1,j1m1,k1m1)
    ind2 = indc(i2+1,j2+1,k2+1)
    do n=ind1,ind2
       m=n-n0c
       d(n,1)=0.
       d(n,2)=0.
       d(n,3)=0.
       d(n,4)=0.
       d(n,5)=0.
       dfxx(m)=0.
       dfxy(m)=0.
       dfxz(m)=0.
       dfex(m)=0.
       dfyy(m)=0.
       dfyz(m)=0.
       dfey(m)=0.
       dfzz(m)=0.
       dfez(m)=0.
       coefe(1,m)=0.
       coefe(2,m)=0.
       coefe(3,m)=0.
       rv(m)=0.
    enddo

!-----rayon spectral visqueux et coef diagonal------------------------------

    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             rv(m)=gam*(mu(n)/pr+mut(n)/prt)/v(n,1)
!         rv(m)=max(gam*(mu(n)/pr+mut(n)/prt)/v(n,1), &
!                   4./3.*(mu(n)+mut(n))/v(n,1))
             coefdiag(m)=vol(n)/dt(n)
          enddo
       enddo
    enddo

!-----remplissage du coefficient diagonal par direction--------------------------

    kdir=1
    ninc=nci

    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1,j,k)
          ind2 = indc(i2,j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
             cc=0.5*(cson(n)+cson(n-ninc))
             coefe(kdir,m)=0.5*(abs(vn) + sqrt(cnds)*cc) &
                  + (rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
          enddo
       enddo
    enddo

    kdir=2
    ninc=ncj

    do k=k1,k2m1
       do j=j1,j2
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
             cc=0.5*(cson(n)+cson(n-ninc))
             coefe(kdir,m)=0.5*(abs(vn) + sqrt(cnds)*cc) &
                  +(rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
          enddo
       enddo
    enddo

    kdir=3
    ninc=nck

    do k=k1,k2
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             vn=uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3)
             cc=0.5*(cson(n)+cson(n-ninc))
             coefe(kdir,m)=0.5*(abs(vn) + sqrt(cnds)*cc) &
                  +(rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
          enddo
       enddo
    enddo

    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             coefdiag(m)=coefdiag(m) + coefe(1,m) + coefe(1,m+nci) &
                  + coefe(2,m) + coefe(2,m+ncj) &
                  + coefe(3,m) + coefe(3,m+nck)
          enddo
       enddo
    enddo

!------calcul instationnaire avec dts-----------------------------------

    if(kfmg.eq.3) then
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                coefdiag(m)=coefdiag(m) + fact*vol(n)/dtpas
             enddo
          enddo
       enddo
    endif

!*************************************************************************
!    boucle sur les sous-iterations
!*************************************************************************
!
    do ls=1,lmx

!-----residu explicite------------------------------------------

       if(ityprk.eq.0) then
          do k=k1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   d2w1(m)=-u(n,1)
                   d2w2(m)=-u(n,2)
                   d2w3(m)=-u(n,3)
                   d2w4(m)=-u(n,4)
                   d2w5(m)=-u(n,5)
                enddo
             enddo
          enddo
       else
          do k=k1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   d2w1(m)=-u(n,1)-ff(n,1)
                   d2w2(m)=-u(n,2)-ff(n,2)
                   d2w3(m)=-u(n,3)-ff(n,3)
                   d2w4(m)=-u(n,4)-ff(n,4)
                   d2w5(m)=-u(n,5)-ff(n,5)
                enddo
             enddo
          enddo
       endif

!------direction i------------------------------------------

       kdir=1
       ninc=nci

       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1,j,k)
             ind2 = indc(i2,j,k)
             do n=ind1,ind2
                m=n-n0c
                tn1=0.5*(d(n,2)+d(n-ninc,2))*sn(m,kdir,1) &
                   +0.5*(d(n,3)+d(n-ninc,3))*sn(m,kdir,2) &
                   +0.5*(d(n,4)+d(n-ninc,4))*sn(m,kdir,3)
                tn2=0.5*(dfxx(m)+dfxx(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfxy(m)+dfxy(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfxz(m)+dfxz(m-ninc))*sn(m,kdir,3)
                tn3=0.5*(dfxy(m)+dfxy(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfyy(m)+dfyy(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfyz(m)+dfyz(m-ninc))*sn(m,kdir,3)
                tn4=0.5*(dfxz(m)+dfxz(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfyz(m)+dfyz(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfzz(m)+dfzz(m-ninc))*sn(m,kdir,3)
                tn5=0.5*(dfex(m)+dfex(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfey(m)+dfey(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfez(m)+dfez(m-ninc))*sn(m,kdir,3)
                d2w1(m)=d2w1(m) + tn1 &
                     + coefe(kdir,m)*d(n-ninc,1) &
                     + coefe(kdir,m+ninc)*d(n+ninc,1)
                d2w2(m)=d2w2(m) + tn2 &
                     + coefe(kdir,m)*d(n-ninc,2) &
                     + coefe(kdir,m+ninc)*d(n+ninc,2)
                d2w3(m)=d2w3(m) + tn3 &
                     + coefe(kdir,m)*d(n-ninc,3) &
                     + coefe(kdir,m+ninc)*d(n+ninc,3)
                d2w4(m)=d2w4(m) + tn4 &
                     + coefe(kdir,m)*d(n-ninc,4) &
                     + coefe(kdir,m+ninc)*d(n+ninc,4)
                d2w5(m)=d2w5(m) + tn5 &
                     + coefe(kdir,m)*d(n-ninc,5) &
                     + coefe(kdir,m+ninc)*d(n+ninc,5)
                d2w1(m-ninc)=d2w1(m-ninc) - tn1
                d2w2(m-ninc)=d2w2(m-ninc) - tn2
                d2w3(m-ninc)=d2w3(m-ninc) - tn3
                d2w4(m-ninc)=d2w4(m-ninc) - tn4
                d2w5(m-ninc)=d2w5(m-ninc) - tn5
             enddo
          enddo
       enddo

!------direction j------------------------------------------

       kdir=2
       ninc=ncj

       do k=k1,k2m1
          do j=j1,j2
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                tn1=0.5*(d(n,2)+d(n-ninc,2))*sn(m,kdir,1) &
                     +0.5*(d(n,3)+d(n-ninc,3))*sn(m,kdir,2) &
                     +0.5*(d(n,4)+d(n-ninc,4))*sn(m,kdir,3)
                tn2=0.5*(dfxx(m)+dfxx(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfxy(m)+dfxy(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfxz(m)+dfxz(m-ninc))*sn(m,kdir,3)
                tn3=0.5*(dfxy(m)+dfxy(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfyy(m)+dfyy(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfyz(m)+dfyz(m-ninc))*sn(m,kdir,3)
                tn4=0.5*(dfxz(m)+dfxz(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfyz(m)+dfyz(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfzz(m)+dfzz(m-ninc))*sn(m,kdir,3)
                tn5=0.5*(dfex(m)+dfex(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfey(m)+dfey(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfez(m)+dfez(m-ninc))*sn(m,kdir,3)
                d2w1(m)=d2w1(m) + tn1 &
                     + coefe(kdir,m)*d(n-ninc,1) &
                     + coefe(kdir,m+ninc)*d(n+ninc,1)
                d2w2(m)=d2w2(m) + tn2 &
                     + coefe(kdir,m)*d(n-ninc,2) &
                     + coefe(kdir,m+ninc)*d(n+ninc,2)
                d2w3(m)=d2w3(m) + tn3 &
                     + coefe(kdir,m)*d(n-ninc,3) &
                     + coefe(kdir,m+ninc)*d(n+ninc,3)
                d2w4(m)=d2w4(m) + tn4 &
                     + coefe(kdir,m)*d(n-ninc,4) &
                     + coefe(kdir,m+ninc)*d(n+ninc,4)
                d2w5(m)=d2w5(m) + tn5 &
                     + coefe(kdir,m)*d(n-ninc,5) &
                     + coefe(kdir,m+ninc)*d(n+ninc,5)
                d2w1(m-ninc)=d2w1(m-ninc) - tn1
                d2w2(m-ninc)=d2w2(m-ninc) - tn2
                d2w3(m-ninc)=d2w3(m-ninc) - tn3
                d2w4(m-ninc)=d2w4(m-ninc) - tn4
                d2w5(m-ninc)=d2w5(m-ninc) - tn5
             enddo
          enddo
       enddo

!------direction k------------------------------------------

       kdir=3
       ninc=nck

       do k=k1,k2
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                tn1=0.5*(d(n,2)+d(n-ninc,2))*sn(m,kdir,1) &
                     +0.5*(d(n,3)+d(n-ninc,3))*sn(m,kdir,2) &
                     +0.5*(d(n,4)+d(n-ninc,4))*sn(m,kdir,3)
                tn2=0.5*(dfxx(m)+dfxx(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfxy(m)+dfxy(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfxz(m)+dfxz(m-ninc))*sn(m,kdir,3)
                tn3=0.5*(dfxy(m)+dfxy(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfyy(m)+dfyy(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfyz(m)+dfyz(m-ninc))*sn(m,kdir,3)
                tn4=0.5*(dfxz(m)+dfxz(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfyz(m)+dfyz(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfzz(m)+dfzz(m-ninc))*sn(m,kdir,3)
                tn5=0.5*(dfex(m)+dfex(m-ninc))*sn(m,kdir,1) &
                     +0.5*(dfey(m)+dfey(m-ninc))*sn(m,kdir,2) &
                     +0.5*(dfez(m)+dfez(m-ninc))*sn(m,kdir,3)
                d2w1(m)=d2w1(m) + tn1 &
                     + coefe(kdir,m)*d(n-ninc,1) &
                     + coefe(kdir,m+ninc)*d(n+ninc,1)
                d2w2(m)=d2w2(m) + tn2 &
                     + coefe(kdir,m)*d(n-ninc,2) &
                     + coefe(kdir,m+ninc)*d(n+ninc,2)
                d2w3(m)=d2w3(m) + tn3 &
                     + coefe(kdir,m)*d(n-ninc,3) &
                     + coefe(kdir,m+ninc)*d(n+ninc,3)
                d2w4(m)=d2w4(m) + tn4 &
                     + coefe(kdir,m)*d(n-ninc,4) &
                     + coefe(kdir,m+ninc)*d(n+ninc,4)
                d2w5(m)=d2w5(m) + tn5 &
                     + coefe(kdir,m)*d(n-ninc,5) &
                     + coefe(kdir,m+ninc)*d(n+ninc,5)
                d2w1(m-ninc)=d2w1(m-ninc) - tn1
                d2w2(m-ninc)=d2w2(m-ninc) - tn2
                d2w3(m-ninc)=d2w3(m-ninc) - tn3
                d2w4(m-ninc)=d2w4(m-ninc) - tn4
                d2w5(m-ninc)=d2w5(m-ninc) - tn5
             enddo
          enddo
       enddo

!*******************************************************************************
!     Calcul de l'increment implicite
!     Actualisation des variables conservatives et des flux
!     Calcul des increments de flux
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                d(n,1)=d2w1(m)/coefdiag(m)
                d(n,2)=d2w2(m)/coefdiag(m)
                d(n,3)=d2w3(m)/coefdiag(m)
                d(n,4)=d2w4(m)/coefdiag(m)
                d(n,5)=d2w5(m)/coefdiag(m)
!
                wi1=v(n,1)+d(n,1)
                wi2=v(n,2)+d(n,2)
                wi3=v(n,3)+d(n,3)
                wi4=v(n,4)+d(n,4)
                wi5=v(n,5)+d(n,5)
                ui=wi2/wi1
                vi=wi3/wi1
                wi=wi4/wi1
!         pression donnee par loi d'etat
                pres=gam1*(wi5-.5*wi1*(ui**2+vi**2+wi**2)-pinfl)
!
                fxx=v(n,2)*(v(n,2)/v(n,1))+ps(n)
                fxy=v(n,3)*(v(n,2)/v(n,1))
                fxz=v(n,4)*(v(n,2)/v(n,1))
                fyy=v(n,3)*(v(n,3)/v(n,1))+ps(n)
                fyz=v(n,4)*(v(n,3)/v(n,1))
                fzz=v(n,4)*(v(n,4)/v(n,1))+ps(n)
                fex=(v(n,5)+ps(n)-pinfl)*v(n,2)/v(n,1)
                fey=(v(n,5)+ps(n)-pinfl)*v(n,3)/v(n,1)
                fez=(v(n,5)+ps(n)-pinfl)*v(n,4)/v(n,1)
!
                dfxx(m)=wi1*ui**2+pres-fxx
                dfxy(m)=wi1*ui*vi-fxy
                dfxz(m)=wi1*ui*wi-fxz
                dfex(m)=ui*(wi5+pres-pinfl)-fex
                dfyy(m)=wi1*vi**2+pres-fyy
                dfyz(m)=wi1*vi*wi-fyz
                dfey(m)=vi*(wi5+pres-pinfl)-fey
                dfzz(m)=wi1*wi**2+pres-fzz
                dfez(m)=wi*(wi5+pres-pinfl)-fez
             enddo
          enddo
       enddo

    enddo  !fin boucle sous-iterations
!
!*************************************************************************
!      avance d'un pas de temps des variables
!*************************************************************************

    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             v(n,1)=v(n,1)+d(n,1)
             v(n,2)=v(n,2)+d(n,2)
             v(n,3)=v(n,3)+d(n,3)
             v(n,4)=v(n,4)+d(n,4)
             v(n,5)=v(n,5)+d(n,5)
          enddo
       enddo
    enddo

    DEALLOCATE(coefe,coefdiag,d2w1,d2w2,d2w3,d2w4,d2w5)

    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine implimf_3d
end module mod_implimf_3d
