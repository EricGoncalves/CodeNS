module mod_implimf_eu
  implicit none
contains
  subroutine implimf_eu( &
       lm,u,dt,v,d,ff, &
       equat,lmx, &
       sn,lgsnlt, &
       vol,dtpas,ityprk, &
       dfxx,dfyy,dfzz,dfxy,dfxz,dfyz,dfex,dfey,dfez,coefdiag, &
       ps,cson)
!
!***********************************************************************
!
!_DA  DATE_C : mars 2002 - Eric Goncalves / Sinumef
!
!     ACT
!_A    Phase implicite sans matrice avec relaxation
!_A    au moyen d'une methode de Jacobi par points.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use schemanum
    implicit none
  integer          ::      i,    i1,  i1m1,    i2,  i2m1
  integer          ::     id,  ind1,  ind2,ityprk,     j
  integer          ::     j1,  j1m1,    j2,  j2m1,    jd
  integer          ::      k,    k1,  k1m1,    k2,  k2m1
  integer          ::     kd,  kdir,lgsnlt,    lm,   lmx
  integer          ::     ls,     m,     n,   n0c,   nci
  integer          ::    ncj,   nck,   nid,  nijd,  ninc
  integer          ::    njd
  double precision ::                   cc,                cnds,      coefdiag(ip00),          cson(ip11),        d(ip11,ip60)
  double precision ::           dfex(ip00),          dfey(ip00),          dfez(ip00),          dfxx(ip00),          dfxy(ip00)
  double precision ::           dfxz(ip00),          dfyy(ip00),          dfyz(ip00),          dfzz(ip00),            dt(ip11)
  double precision ::                dtpas,                fact,                 fex,                 fey,                 fez
  double precision ::        ff(ip11,ip60),                fiex,                fiey,                fiez,                fixx
  double precision ::                 fixy,                fixz,                fiyy,                fiyz,                fizz
  double precision ::                  fxx,                 fxy,                 fxz,                 fyy,                 fyz
  double precision ::                  fzz,                pres,            ps(ip11),sn(lgsnlt,nind,ndir),                 tn1
  double precision ::                  tn2,                 tn3,                 tn4,                 tn5,        u(ip11,ip60)
  double precision ::                   ui,                  uu,        v(ip11,ip60),                  vi,           vol(ip11)
  double precision ::                   vv,                  wi,                 wi1,                 wi2,                 wi3
  double precision ::                  wi4,                 wi5,                  ww
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!



    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: coefe
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE   :: d2w1,d2w2,d2w3,d2w4,d2w5
    ALLOCATE(coefe(ndir,ip00))
    ALLOCATE(d2w1(ip00),d2w2(ip00),d2w3(ip00),d2w4(ip00),d2w5(ip00))

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
!
!-----initalisation--------------------------------
!
    ind1 = indc(i1m1,j1m1,k1m1)
    ind2 = indc(i2+1,j2+1,k2+1)
!!!!$OMP PARALLEL default(SHARED) 
!!!!$OMP DO 
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
    enddo
!!!!$OMP END DO 
!
!------coef diagonal ------------------------------------------------
!
    do k=k1,k2m1
!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             coefdiag(m)=vol(n)/dt(n)
          enddo
       enddo
!!!!$OMP END DO
    enddo
!
!-----remplissage du coefficient diagonal par direction---------------
!
!!!!$OMP SINGLE
    kdir=1
    ninc=nci
!!!!$OMP END SINGLE
!
    do k=k1,k2m1
!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,cnds,uu,vv,ww,cc)
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
             cc=0.5*(cson(n)+cson(n-ninc))
             coefe(kdir,m)=0.5*(abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2) &
                  + ww*sn(m,kdir,3))) + sqrt(cnds)*cc
          enddo
       enddo
!!!!$OMP END DO
    enddo
!
!!!!$OMP SINGLE
    kdir=2
    ninc=ncj
!!!!$OMP END SINGLE
!
    do k=k1,k2m1
!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,cnds,uu,vv,ww,cc)
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
             cc=0.5*(cson(n)+cson(n-ninc))
             coefe(kdir,m)=0.5*(abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2) &
                  + ww*sn(m,kdir,3))) + sqrt(cnds)*cc
          enddo
       enddo
!!!!$OMP END DO
    enddo
!
    do k=k1,k2m1
!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             coefdiag(m)=coefdiag(m) + coefe(1,m) + coefe(1,m+nci) &
                  + coefe(2,m) + coefe(2,m+ncj)
          enddo
       enddo
!!!!$OMP END DO
    enddo
!
!      calcul instationnaire avec pas de temps dual
!
    if(kfmg.eq.3) then
       fact=1.5
       do k=k1,k2m1
!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                coefdiag(m)=coefdiag(m) + fact*vol(n)/dt1min
             enddo
          enddo
!!!!$OMP END DO
       enddo
    endif
!
    if(equat(3:4).eq.'3d') then
!
       kdir=3
       ninc=nck
!
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
                cc=0.5*(cson(n)+cson(n-ninc))
                coefe(kdir,m)=0.5*(abs(uu*sn(m,kdir,1)+vv*sn(m,kdir,2) &
                     + ww*sn(m,kdir,3))) + sqrt(cnds)*cc
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                coefdiag(m)=coefdiag(m) + coefe(3,m) + coefe(3,m+nck)
             enddo
          enddo
       enddo
!
    endif
!
!*************************************************************************
!c    boucle sur les sous-iterations
!*************************************************************************
!
    do ls=1,lmx
!
!-----residu explicite------------------------------------------
!
       if(ityprk.eq.0) then
          do k=k1,k2m1
!!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
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
!!!!!$OMP END DO
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
!
!------direction i------------------------------------------
!
       kdir=1
       ninc=nci
!
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
!
!------direction j------------------------------------------
!
       kdir=2
       ninc=ncj
!
       do k=k1,k2m1
!!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,tn1,tn2,tn3,tn5)
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
!!!!!$OMP END DO
       enddo
!
!------direction k------------------------------------------
!
       if(equat(3:4).eq.'3d') then
          kdir=3
          ninc=nck
!
          do k=k1,k2
!!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,tn1,tn2,tn3,tn5)
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
!!!!!$OMP END DO
          enddo
       endif
!
!*******************************************************************************
!
!c    calcul de l'increment implicite
!c    actualisation des variables conservatives et des flux
!c    calcul des increments de flux
!
       do k=k1,k2m1
!!!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,wi1,wi2,wi3,wi4,wi5,ui,vi,wi,pres,&
!!!!!$OMP fixx,fixy,fixz,fiyy,fiyz,fizz,fiex,fiey,fiez,&
!!!!!$OMP fxx,fxy,fxz,fyy,fyz,fzz,fex,fey,fez)
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
                fixx=wi1*ui**2+pres
                fixy=wi1*ui*vi
                fixz=wi1*ui*wi
                fiyy=wi1*vi**2+pres
                fiyz=wi1*vi*wi
                fizz=wi1*wi**2+pres
                fiex=ui*(wi5+pres-pinfl)
                fiey=vi*(wi5+pres-pinfl)
                fiez=wi*(wi5+pres-pinfl)
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
                dfxx(m)=fixx-fxx
                dfxy(m)=fixy-fxy
                dfxz(m)=fixz-fxz
                dfex(m)=fiex-fex
                dfyy(m)=fiyy-fyy
                dfyz(m)=fiyz-fyz
                dfey(m)=fiey-fey
                dfzz(m)=fizz-fzz
                dfez(m)=fiez-fez
             enddo
          enddo
!!!!!$OMP END DO
       enddo
!
    enddo  !fin boucle sous-iterations
!
!*************************************************************************
!      avance d'un pas de temps des variables
!*************************************************************************
!
!
    do k=k1,k2m1
!!!!!$OMP DO PRIVATE(j,n,ind1,ind2)
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             v(n,:)=v(n,:)+d(n,:)
          enddo
       enddo
!!!!!$OMP END DO
    enddo
!!!!$OMP END PARALLEL

    DEALLOCATE(coefe,d2w1,d2w2,d2w3,d2w4,d2w5)

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
  end subroutine implimf_eu
end module mod_implimf_eu
