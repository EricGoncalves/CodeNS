module mod_implimf_prcd2_3d
  implicit none
contains
  subroutine implimf_prcd2_3d( &
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
!_DA  DATE_C : janvier 2007 - Eric Goncalves / LEGI
!
!     ACT
!_A    Phase implicite sans matrice avec relaxation
!_A    au moyen d'une methode de Jacobi par points.
!_A    Preconditionnement basse vitesse de Turkel (P,u,e)
!_A    Pour calculs paralleles en geometrie 3D.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use schemanum
    use definition 
    implicit none
    integer          ::      i,    i1,  i1m1,    i2,  i2m1
    integer          ::     id,  ind1,  ind2,ityprk,     j
    integer          ::     j1,  j1m1,    j2,  j2m1,    jd
    integer          ::      k,    k1,  k1m1,    k2,  k2m1
    integer          ::     kd,  kdir,lgsnlt,    lm,   lmx
    integer          ::     ls,     m,     n,   n0c,   nci
    integer          ::    ncj,   nck,   nid,  nijd,  ninc
    integer          ::    njd
    double precision ::     a2, beta2,  cnds, cndsi, cndsj
    double precision ::  coefa, cson(ip11), d(ip11,ip60), dfex(ip00),dfey(ip00)
    double precision ::  dfez(ip00), dfxx(ip00), dfxy(ip00), dfxz(ip00),dfyy(ip00)
    double precision ::  dfyz(ip00), dfzz(ip00), dt(ip11), dtpas, cndsk 
    double precision ::   dw11,  dw12,  dw13,  dw14,  dw15  
    double precision ::   dw21,  dw22,  dw23,  dw24,  dw25
    double precision ::   fact,fex,fey,sn(lgsnlt,nind,ndir), u(ip11,ip60)
    double precision ::    fez, ff(ip11,ip60), fxx, fxy, fxz
    double precision ::    fyy,   fyz,   fzz,    gd,    ge
    double precision ::    get, mu(ip12), mut(ip12), precon,  pres
    double precision ::     ps(ip11), q2,  qinf,  rhoe, rv(ip00)
    double precision ::    ti1,   ti2,   ti3,   ti4,   ti5
    double precision ::    tj1,   tj2,   tj3,   tj4,   tj5
    double precision ::    tk1,   tk2,   tk3,   tk4,   tk5
    double precision ::     ui,    uu, v(ip11,ip60),    vi,    vn
    double precision ::    vol(ip11), vv, wi,   wi1,   wi2
    double precision ::    wi3,   wi4,   wi5,    ww
    double precision,allocatable :: coefb(:),coefdiag(:),coefe(:,:),coefv(:,:)
    double precision,allocatable :: d2w1(:),d2w2(:),d2w3(:),d2w4(:),d2w5(:)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat

    ALLOCATE(coefe(ndir,ip00),coefv(ndir,ip00))
    ALLOCATE(coefdiag(ip00),coefb(ip00), &
            d2w1(ip00),d2w2(ip00),d2w3(ip00),d2w4(ip00),d2w5(ip00))
!
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
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
      fact=1.5
!
!-----initalisation--------------------------------
!
      ind1 = indc(i1m1,j1m1,k1m1)
      ind2 = indc(i2+1,j2+1,k2+1)
!!!$OMP PARALLEL 
!!!$OMP DO 
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
       dfyy(m)=0.
       dfyz(m)=0.
       dfzz(m)=0.
       dfex(m)=0.
       dfey(m)=0.
       dfez(m)=0.
       coefe(1,m)=0.
       coefe(2,m)=0.
       coefe(3,m)=0.
       coefv(1,m)=0.
       coefv(2,m)=0.
       coefv(3,m)=0.
       rv(m)=0.
      enddo
!!!$OMP END DO 
!
!------rayon spectral visqueux et coef diagonal------------------------------
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2)
      do k=k1,k2m1
       do j=j1,j2m1
        ind1 = indc(i1  ,j,k)
        ind2 = indc(i2m1,j,k)
        do n=ind1,ind2
         m=n-n0c
         rv(m)=max(gam*(mu(n)/pr+mut(n)/prt)/v(n,1), &
                   4./3.*(mu(n)+mut(n))/v(n,1))
         coefdiag(m)=vol(n)/dt(n)
         coefb(m)=0.
        enddo
       enddo
      enddo
!!!$OMP END DO
!
!-----remplissage du coefficient diagonal par direction------------------
!
!!!$OMP SINGLE
       kdir=1
       ninc=nci
!!!$OMP END SINGLE
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2,cnds,uu,vv,ww,vn,a2,beta2)
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
          a2=0.25*(cson(n)+cson(n-ninc))**2
          beta2=min(max((uu**2+vv**2+ww**2)/a2,cte*qinf**2/a2),1.)
          coefe(kdir,m)=0.25*((1.+beta2)*abs(vn) &
                     +    sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          coefv(kdir,m)=(rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
         enddo
        enddo
       enddo
!!!$OMP END DO
!
!!!$OMP SINGLE
       kdir=2
       ninc=ncj
!!!$OMP END SINGLE
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2,cnds,uu,vv,ww,vn,a2,beta2)
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
          a2=0.25*(cson(n)+cson(n-ninc))**2
          beta2=min(max((uu**2+vv**2+ww**2)/a2,cte*qinf**2/a2),1.)
          coefe(kdir,m)=0.25*((1.+beta2)*abs(vn) &
                     +    sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          coefv(kdir,m)=(rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
         enddo
        enddo
       enddo
!!!$OMP END DO
!
!!!$OMP SINGLE
       kdir=3
       ninc=nck
!!!$OMP END SINGLE
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2,cnds,uu,vv,ww,vn,a2,beta2)
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
          a2=0.25*(cson(n)+cson(n-ninc))**2
          beta2=min(max((uu**2+vv**2+ww**2)/a2,cte*qinf**2/a2),1.)
          coefe(kdir,m)=0.25*((1.+beta2)*abs(vn) &
                     +    sqrt(((1.-beta2)*vn)**2+4.*beta2*cnds*a2))
          coefv(kdir,m)=(rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
         enddo
        enddo
       enddo
!!!$OMP END DO
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2)
       do k=k1,k2m1
        do j=j1,j2m1
         ind1 = indc(i1  ,j,k)
         ind2 = indc(i2m1,j,k)
         do n=ind1,ind2
          m=n-n0c
          coefb(m)=coefv(1,m) + coefv(1,m+nci) &
                 + coefv(2,m) + coefv(2,m+ncj) &
                 + coefv(3,m) + coefv(3,m+nck)
          coefdiag(m)=coefdiag(m) + coefb(m) &
                     +coefe(1,m) + coefe(1,m+nci) &
                     +coefe(2,m) + coefe(2,m+ncj) &
                     +coefe(3,m) + coefe(3,m+nck)
         enddo
        enddo
       enddo
!!!$OMP END DO
!
!-----calcul instationnaire avec dts-------------------------------------
!
       if(kfmg.eq.3) then
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2)
        do k=k1,k2m1
         do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
           m=n-n0c
           coefdiag(m)=coefdiag(m) + fact*vol(n)/dt1min
           coefb(m)   =coefb(m)    + fact*vol(n)/dt1min
          enddo
         enddo
        enddo
!!!$OMP END DO
       endif
!
!*************************************************************************
!    boucle sur les sous-iterations
!*************************************************************************
!
      do ls=1,lmx
!
!-----residu explicite------------------------------------------
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2)
        do k=k1,k2m1
         do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
!DEC$ IVDEP
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
!!!$OMP END DO
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2,ti1,ti2,ti3,ti4,ti5,tj1,tj2,tj3,tj4,tj5,tk1,tk2,tk3,tk4,tk5,cndsi,cndsj,cndsk, &
!!!$OMP uu,vv,ww,q2,a2,beta2,get,ge,coefa,gd,dw11,dw12,dw13,dw14,dw15,dw21,dw22,dw23,dw24,dw25,precon)
       do k=k1,k2m1
        do j=j1,j2m1
         ind1 = indc(i1  ,j,k)
         ind2 = indc(i2m1,j,k)
         do n=ind1,ind2
          m=n-n0c
          ti1=(d(n,2)+d(n-nci,2))*sn(m,1,1) &
             +(d(n,3)+d(n-nci,3))*sn(m,1,2) &
             +(d(n,4)+d(n-nci,4))*sn(m,1,3) &
             -(d(n,2)+d(n+nci,2))*sn(m+nci,1,1) &
             -(d(n,3)+d(n+nci,3))*sn(m+nci,1,2) &
             -(d(n,4)+d(n+nci,4))*sn(m+nci,1,3)
          ti2=(dfxx(m)+dfxx(m-nci))*sn(m,1,1) &
             +(dfxy(m)+dfxy(m-nci))*sn(m,1,2) &
             +(dfxz(m)+dfxz(m-nci))*sn(m,1,3) &
             -(dfxx(m)+dfxx(m+nci))*sn(m+nci,1,1) &
             -(dfxy(m)+dfxy(m+nci))*sn(m+nci,1,2) &
             -(dfxz(m)+dfxz(m+nci))*sn(m+nci,1,3)
          ti3=(dfxy(m)+dfxy(m-nci))*sn(m,1,1) &
             +(dfyy(m)+dfyy(m-nci))*sn(m,1,2) &
             +(dfyz(m)+dfyz(m-nci))*sn(m,1,3) &
             -(dfxy(m)+dfxy(m+nci))*sn(m+nci,1,1) &
             -(dfyy(m)+dfyy(m+nci))*sn(m+nci,1,2) &
             -(dfyz(m)+dfyz(m+nci))*sn(m+nci,1,3)
          ti4=(dfxz(m)+dfxz(m-nci))*sn(m,1,1) &
             +(dfyz(m)+dfyz(m-nci))*sn(m,1,2) &
             +(dfzz(m)+dfzz(m-nci))*sn(m,1,3) &
             -(dfxz(m)+dfxz(m+nci))*sn(m+nci,1,1) &
             -(dfyz(m)+dfyz(m+nci))*sn(m+nci,1,2) &
             -(dfzz(m)+dfzz(m+nci))*sn(m+nci,1,3)
          ti5=(dfex(m)+dfex(m-nci))*sn(m,1,1) &
             +(dfey(m)+dfey(m-nci))*sn(m,1,2) &
             +(dfez(m)+dfez(m-nci))*sn(m,1,3) &
             -(dfex(m)+dfex(m+nci))*sn(m+nci,1,1) &
             -(dfey(m)+dfey(m+nci))*sn(m+nci,1,2) &
             -(dfez(m)+dfez(m+nci))*sn(m+nci,1,3)
!
          tj1=(d(n,2)+d(n-ncj,2))*sn(m,2,1) &
             +(d(n,3)+d(n-ncj,3))*sn(m,2,2) &
             +(d(n,4)+d(n-ncj,4))*sn(m,2,3) &
             -(d(n,2)+d(n+ncj,2))*sn(m+ncj,2,1) &
             -(d(n,3)+d(n+ncj,3))*sn(m+ncj,2,2) &
             -(d(n,4)+d(n+ncj,4))*sn(m+ncj,2,3)
          tj2=(dfxx(m)+dfxx(m-ncj))*sn(m,2,1) &
             +(dfxy(m)+dfxy(m-ncj))*sn(m,2,2) &
             +(dfxz(m)+dfxz(m-ncj))*sn(m,2,3) &
             -(dfxx(m)+dfxx(m+ncj))*sn(m+ncj,2,1) &
             -(dfxy(m)+dfxy(m+ncj))*sn(m+ncj,2,2) &
             -(dfxz(m)+dfxz(m+ncj))*sn(m+ncj,2,3)
          tj3=(dfxy(m)+dfxy(m-ncj))*sn(m,2,1) &
             +(dfyy(m)+dfyy(m-ncj))*sn(m,2,2) &
             +(dfyz(m)+dfyz(m-ncj))*sn(m,2,3) &
             -(dfxy(m)+dfxy(m+ncj))*sn(m+ncj,2,1) &
             -(dfyy(m)+dfyy(m+ncj))*sn(m+ncj,2,2) &
             -(dfyz(m)+dfyz(m+ncj))*sn(m+ncj,2,3)
          tj4=(dfxz(m)+dfxz(m-ncj))*sn(m,2,1) &
             +(dfyz(m)+dfyz(m-ncj))*sn(m,2,2) &
             +(dfzz(m)+dfzz(m-ncj))*sn(m,2,3) &
             -(dfxz(m)+dfxz(m+ncj))*sn(m+ncj,2,1) &
             -(dfyz(m)+dfyz(m+ncj))*sn(m+ncj,2,2) &
             -(dfzz(m)+dfzz(m+ncj))*sn(m+ncj,2,3)
          tj5=(dfex(m)+dfex(m-ncj))*sn(m,2,1) &
             +(dfey(m)+dfey(m-ncj))*sn(m,2,2) &
             +(dfez(m)+dfez(m-ncj))*sn(m,2,3) &
             -(dfex(m)+dfex(m+ncj))*sn(m+ncj,2,1) &
             -(dfey(m)+dfey(m+ncj))*sn(m+ncj,2,2) &
             -(dfez(m)+dfez(m+ncj))*sn(m+ncj,2,3)
!
          tk1=(d(n,2)+d(n-nck,2))*sn(m,3,1) &
             +(d(n,3)+d(n-nck,3))*sn(m,3,2) &
             +(d(n,4)+d(n-nck,4))*sn(m,3,3) &
             -(d(n,2)+d(n+nck,2))*sn(m+nck,3,1) &
             -(d(n,3)+d(n+nck,3))*sn(m+nck,3,2) &
             -(d(n,4)+d(n+nck,4))*sn(m+nck,3,3)
          tk2=(dfxx(m)+dfxx(m-nck))*sn(m,3,1) &
             +(dfxy(m)+dfxy(m-nck))*sn(m,3,2) &
             +(dfxz(m)+dfxz(m-nck))*sn(m,3,3) &
             -(dfxx(m)+dfxx(m+nck))*sn(m+nck,3,1) &
             -(dfxy(m)+dfxy(m+nck))*sn(m+nck,3,2) &
             -(dfxz(m)+dfxz(m+nck))*sn(m+nck,3,3)
          tk3=(dfxy(m)+dfxy(m-nck))*sn(m,3,1) &
             +(dfyy(m)+dfyy(m-nck))*sn(m,3,2) &
             +(dfyz(m)+dfyz(m-nck))*sn(m,3,3) &
             -(dfxy(m)+dfxy(m+nck))*sn(m+nck,3,1) &
             -(dfyy(m)+dfyy(m+nck))*sn(m+nck,3,2) &
             -(dfyz(m)+dfyz(m+nck))*sn(m+nck,3,3)
          tk4=(dfxz(m)+dfxz(m-nck))*sn(m,3,1) &
             +(dfyz(m)+dfyz(m-nck))*sn(m,3,2) &
             +(dfzz(m)+dfzz(m-nck))*sn(m,3,3) &
             -(dfxz(m)+dfxz(m+nck))*sn(m+nck,3,1) &
             -(dfyz(m)+dfyz(m+nck))*sn(m+nck,3,2) &
             -(dfzz(m)+dfzz(m+nck))*sn(m+nck,3,3)
          tk5=(dfex(m)+dfex(m-nck))*sn(m,3,1) &
             +(dfey(m)+dfey(m-nck))*sn(m,3,2) &
             +(dfez(m)+dfez(m-nck))*sn(m,3,3) &
             -(dfex(m)+dfex(m+nck))*sn(m+nck,3,1) &
             -(dfey(m)+dfey(m+nck))*sn(m+nck,3,2) &
             -(dfez(m)+dfez(m+nck))*sn(m+nck,3,3)
!
          cndsi=sqrt(sn(m,1,1)*sn(m,1,1)+ &
                     sn(m,1,2)*sn(m,1,2)+ &
                     sn(m,1,3)*sn(m,1,3))
          cndsj=sqrt(sn(m,2,1)*sn(m,2,1)+ &
                     sn(m,2,2)*sn(m,2,2)+ &
                     sn(m,2,3)*sn(m,2,3))
          cndsk=sqrt(sn(m,3,1)*sn(m,3,1)+ &
                     sn(m,3,2)*sn(m,3,2)+ &
                     sn(m,3,3)*sn(m,3,3))
          uu=(v(n,2)*sn(m,1,1)+v(n,3)*sn(m,1,2)+v(n,4)*sn(m,1,3))/(v(n,1)*cndsi)
          vv=(v(n,2)*sn(m,2,1)+v(n,3)*sn(m,2,2)+v(n,4)*sn(m,2,3))/(v(n,1)*cndsj)
          ww=(v(n,2)*sn(m,3,1)+v(n,3)*sn(m,3,2)+v(n,4)*sn(m,3,3))/(v(n,1)*cndsk)
          q2=uu**2+vv**2+ww**2
          a2=cson(n)**2
          beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
          get=v(n,5)/v(n,1)             !energie totale
          ge=get-0.5*q2                 !energie interne
          coefa=coefdiag(m)-coefb(m)
          gd=(beta2-1.)/(ge*(coefa+beta2*coefb(m)))
!
          dw11=d2w1(m)+0.5*(ti1+tj1+tk1)+coefv(1,m    )*d(n-nci,1) &
                                   + coefv(1,m+nci)*d(n+nci,1) &
                                   + coefv(2,m    )*d(n-ncj,1) &
                                   + coefv(2,m+ncj)*d(n+ncj,1) &
                                   + coefv(3,m    )*d(n-nck,1) &
                                   + coefv(3,m+nck)*d(n+nck,1)
          dw12=d2w2(m)+0.5*(ti2+tj2+tk2)+coefv(1,m    )*d(n-nci,2) &
                                   + coefv(1,m+nci)*d(n+nci,2) &
                                   + coefv(2,m    )*d(n-ncj,2) &
                                   + coefv(2,m+ncj)*d(n+ncj,2) &
                                   + coefv(3,m    )*d(n-nck,2) &
                                   + coefv(3,m+nck)*d(n+nck,2)
          dw13=d2w3(m)+0.5*(ti3+tj3+tk3)+coefv(1,m    )*d(n-nci,3) &
                                   + coefv(1,m+nci)*d(n+nci,3) &
                                   + coefv(2,m    )*d(n-ncj,3) &
                                   + coefv(2,m+ncj)*d(n+ncj,3) &
                                   + coefv(3,m    )*d(n-nck,3) &
                                   + coefv(3,m+nck)*d(n+nck,3)
          dw14=d2w4(m)+0.5*(ti4+tj4+tk4)+coefv(1,m    )*d(n-nci,4) &
                                   + coefv(1,m+nci)*d(n+nci,4) &
                                   + coefv(2,m    )*d(n-ncj,4) &
                                   + coefv(2,m+ncj)*d(n+ncj,4) &
                                   + coefv(3,m    )*d(n-nck,4) &
                                   + coefv(3,m+nck)*d(n+nck,4)
          dw15=d2w5(m)+0.5*(ti5+tj5+tk5)+coefv(1,m    )*d(n-nci,5) &
                                   + coefv(1,m+nci)*d(n+nci,5) &
                                   + coefv(2,m    )*d(n-ncj,5) &
                                   + coefv(2,m+ncj)*d(n+ncj,5) &
                                   + coefv(3,m    )*d(n-nck,5) &
                                   + coefv(3,m+nck)*d(n+nck,5)
          dw21=coefe(1,m    )*d(n-nci,1) &
              +coefe(1,m+nci)*d(n+nci,1) &
              +coefe(2,m    )*d(n-ncj,1) &
              +coefe(2,m+ncj)*d(n+ncj,1) &
              +coefe(3,m    )*d(n-nck,1) &
              +coefe(3,m+nck)*d(n+nck,1)
          dw22=coefe(1,m    )*d(n-nci,2) &
              +coefe(1,m+nci)*d(n+nci,2) &
              +coefe(2,m    )*d(n-ncj,2) &
              +coefe(2,m+ncj)*d(n+ncj,2) &
              +coefe(3,m    )*d(n-nck,2) &
              +coefe(3,m+nck)*d(n+nck,2)
          dw23=coefe(1,m    )*d(n-nci,3) &
              +coefe(1,m+nci)*d(n+nci,3) &
              +coefe(2,m    )*d(n-ncj,3) &
              +coefe(2,m+ncj)*d(n+ncj,3) &
              +coefe(3,m    )*d(n-nck,3) &
              +coefe(3,m+nck)*d(n+nck,3)
          dw24=coefe(1,m    )*d(n-nci,4) &
              +coefe(1,m+nci)*d(n+nci,4) &
              +coefe(2,m    )*d(n-ncj,4) &
              +coefe(2,m+ncj)*d(n+ncj,4) &
              +coefe(3,m    )*d(n-nck,4) &
              +coefe(3,m+nck)*d(n+nck,4)
          dw25=coefe(1,m    )*d(n-nci,5) &
              +coefe(1,m+nci)*d(n+nci,5) &
              +coefe(2,m    )*d(n-ncj,5) &
              +coefe(2,m+ncj)*d(n+ncj,5) &
              +coefe(3,m    )*d(n-nck,5) &
              +coefe(3,m+nck)*d(n+nck,5)
!
          precon=gd*(0.5*q2*(coefa*dw11-coefb(m)*dw21) &
                        -uu*(coefa*dw12-coefb(m)*dw22) &
                        -vv*(coefa*dw13-coefb(m)*dw23) &
                        -ww*(coefa*dw14-coefb(m)*dw24) &
                           +(coefa*dw15-coefb(m)*dw25))
          d2w1(m)=dw11 + dw21 +    precon
          d2w2(m)=dw12 + dw22 + uu*precon
          d2w3(m)=dw13 + dw23 + vv*precon
          d2w4(m)=dw14 + dw24 + ww*precon
          d2w5(m)=dw15 + dw25 +get*precon
         enddo
        enddo
       enddo
!!!$OMP END DO
!
!*************************************************************************
!    Calcul de l'increment implicite
!    Actualisation des variables conservatives et des flux
!    Calcul des increments de flux
!
!!!$OMP DO PRIVATE(k,j,n,m,ind1,ind2,wi1,wi2,wi3,wi4,wi5,ui,vi,wi,pres,fxx,fxy,fxz,fyy,fyz,fzz,fex,fey,fez)
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
          pres=gam1*(wi5-0.5*wi1*(ui**2+vi**2+wi**2)-pinfl)
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
          dfxy(m)=wi1*ui*vi     -fxy
          dfxz(m)=wi1*ui*wi     -fxz
          dfyy(m)=wi1*vi**2+pres-fyy
          dfyz(m)=wi1*vi*wi     -fyz
          dfzz(m)=wi1*wi**2+pres-fzz
          dfex(m)=ui*(wi5+pres-pinfl) -fex
          dfey(m)=vi*(wi5+pres-pinfl) -fey
          dfez(m)=wi*(wi5+pres-pinfl) -fez
         enddo
        enddo
       enddo
!!!$OMP END DO
!
      enddo  !fin boucle sous-iterations
!
!*************************************************************************
!      avance d'un pas de temps des variables
!*************************************************************************
!
!!!$OMP DO PRIVATE(k,j,n,ind1,ind2)
       do k=k1,k2m1
        do j=j1,j2m1
         ind1 = indc(i1  ,j,k)
         ind2 = indc(i2m1,j,k)
         do n=ind1,ind2
          v(n,1)=v(n,1)+d(n,1)
          v(n,2)=v(n,2)+d(n,2)
          v(n,3)=v(n,3)+d(n,3)
          v(n,4)=v(n,4)+d(n,4)
          v(n,5)=v(n,5)+d(n,5)
        enddo
       enddo
      enddo
!!!$OMP END DO
!!!$OMP END PARALLEL
!
    DEALLOCATE(coefe,coefv,coefdiag,coefb,d2w1,d2w2,d2w3,d2w4,d2w5)

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
  end subroutine implimf_prcd2_3d
end module mod_implimf_prcd2_3d
