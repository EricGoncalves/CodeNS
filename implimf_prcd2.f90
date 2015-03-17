module mod_implimf_prcd2
  implicit none
contains
  subroutine implimf_prcd2( &
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
    integer          ::     id,   inc,  ind1,  ind2,  indc
    integer          :: ityprk,     j,    j1,  j1m1,    j2
    integer          ::   j2m1,    jd,     k,    k1,  k1m1
    integer          ::     k2,  k2m1,    kd,  kdir,lgsnlt
    integer          ::     lm,   lmx,    ls,     m,     n
    integer          ::    n0c,   nci,   ncj,   nck,   nid
    integer          ::   nijd,  ninc,   njd
    double precision ::     a2, beta2,  cnds, cndsi, cndsj
    double precision ::  coefa,  cson,     d,  dfex,  dfey
    double precision ::   dfez,  dfxx,  dfxy,  dfxz,  dfyy
    double precision ::   dfyz,  dfzz,    dt, dtpas,  dw11
    double precision ::   dw12,  dw13,  dw15,  dw21,  dw22
    double precision ::   dw23,  dw25,  fact,   fex,   fey
    double precision ::    fez,    ff,   fxx,   fxy,   fxz
    double precision ::    fyy,   fyz,   fzz,    gd,    ge
    double precision ::    get,    mu,   mut,precon,  pres
    double precision ::     ps,    q2,  qinf,  rhoe,    rv
    double precision ::     sn,   ti1,   ti2,   ti3,   ti5
    double precision ::    tj1,   tj2,   tj3,   tj5,     u
    double precision ::     ui,    uu,     v,    vi,    vn
    double precision ::    vol,    vv,    wi,   wi1,   wi2
    double precision ::    wi3,   wi4,   wi5,    ww
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
    dimension v(ip11,ip60),u(ip11,ip60),d(ip11,ip60),ff(ip11,ip60)
    dimension vol(ip11),dt(ip11),ps(ip11),cson(ip11)
    dimension mu(ip12),mut(ip12)
    dimension sn(lgsnlt,nind,ndir)
    dimension dfxx(ip00),dfxy(ip00),dfxz(ip00),dfex(ip00),rv(ip00), &
         dfyy(ip00),dfyz(ip00),dfey(ip00),dfzz(ip00),dfez(ip00)
!
    indc(i,j,k)=n0c+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    inc(id,jd,kd)=id+jd*nid+kd*nijd

    DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: coefe,coefv
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE   :: coefdiag,coefb, &
         d2w1,d2w2,d2w3,d2w4,d2w5
    ALLOCATE(coefe(ndir,ip00),coefv(ndir,ip00))
    ALLOCATE(coefdiag(ip00),coefb(ip00), &
         d2w1(ip00),d2w2(ip00),d2w3(ip00),d2w4(ip00),d2w5(ip00))

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
    do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
!         rv(m)=(gam/v(n,1))*(mu(n)/pr+mut(n)/prt)
             rv(m)=max(gam*(mu(n)/pr+mut(n)/prt)/v(n,1), &
                  4./3.*(mu(n)+mut(n))/v(n,1))
             coefdiag(m)=vol(n)/dt(n)
             coefb(m)=0.
          enddo
       enddo
!!!$OMP END DO
    enddo
!
!-----remplissage du coefficient diagonal par direction------------------
!
!!!$OMP SINGLE
    kdir=1
    ninc=nci
!!!$OMP END SINGLE
!
    do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,cnds,uu,vv,ww,vn,a2,beta2)
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
!!!$OMP END DO
    enddo
!
!!!$OMP SINGLE
    kdir=2
    ninc=ncj
!!!$OMP END SINGLE
!
    do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,cnds,uu,vv,ww,vn,a2,beta2)
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
!!!$OMP END DO
    enddo
!
    do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             coefb(m)=coefv(1,m) + coefv(1,m+nci) &
                  + coefv(2,m) + coefv(2,m+ncj)
             coefdiag(m)=coefdiag(m) + coefb(m) &
                  +coefe(1,m) + coefe(1,m+nci) &
                  +coefe(2,m) + coefe(2,m+ncj)
          enddo
       enddo
!!!$OMP END DO
    enddo
!
!-----calcul instationnaire avec dts-------------------------------------
!
    if(kfmg.eq.3) then
       fact=1.5
       do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                coefdiag(m)=coefdiag(m) + fact*vol(n)/dt1min
                coefb(m)   =coefb(m)    + fact*vol(n)/dt1min
             enddo
          enddo
!!!$OMP END DO
       enddo
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
       do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2)
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
!!!$OMP END DO
       enddo
!
       do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,ti1,ti2,ti3,ti5,tj1,tj2,tj3,tj5,cndsi,cndsj,uu, &
!!!$OMP vv,q2,a2,beta2,get,ge,coefa,gd,dw11,dw12,dw13,dw15,dw21,dw22,dw23,dw25,precon)
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                ti1=(d(n,2)+d(n-nci,2))*sn(m,1,1) &
                     +(d(n,3)+d(n-nci,3))*sn(m,1,2) &
                     -(d(n,2)+d(n+nci,2))*sn(m+nci,1,1) &
                     -(d(n,3)+d(n+nci,3))*sn(m+nci,1,2)
                ti2=(dfxx(m)+dfxx(m-nci))*sn(m,1,1) &
                     +(dfxy(m)+dfxy(m-nci))*sn(m,1,2) &
                     -(dfxx(m)+dfxx(m+nci))*sn(m+nci,1,1) &
                     -(dfxy(m)+dfxy(m+nci))*sn(m+nci,1,2)
                ti3=(dfxy(m)+dfxy(m-nci))*sn(m,1,1) &
                     +(dfyy(m)+dfyy(m-nci))*sn(m,1,2) &
                     -(dfxy(m)+dfxy(m+nci))*sn(m+nci,1,1) &
                     -(dfyy(m)+dfyy(m+nci))*sn(m+nci,1,2)
                ti5=(dfex(m)+dfex(m-nci))*sn(m,1,1) &
                     +(dfey(m)+dfey(m-nci))*sn(m,1,2) &
                     -(dfex(m)+dfex(m+nci))*sn(m+nci,1,1) &
                     -(dfey(m)+dfey(m+nci))*sn(m+nci,1,2)
!
                tj1=(d(n,2)+d(n-ncj,2))*sn(m,2,1) &
                     +(d(n,3)+d(n-ncj,3))*sn(m,2,2) &
                     -(d(n,2)+d(n+ncj,2))*sn(m+ncj,2,1) &
                     -(d(n,3)+d(n+ncj,3))*sn(m+ncj,2,2)
                tj2=(dfxx(m)+dfxx(m-ncj))*sn(m,2,1) &
                     +(dfxy(m)+dfxy(m-ncj))*sn(m,2,2) &
                     -(dfxx(m)+dfxx(m+ncj))*sn(m+ncj,2,1) &
                     -(dfxy(m)+dfxy(m+ncj))*sn(m+ncj,2,2)
                tj3=(dfxy(m)+dfxy(m-ncj))*sn(m,2,1) &
                     +(dfyy(m)+dfyy(m-ncj))*sn(m,2,2) &
                     -(dfxy(m)+dfxy(m+ncj))*sn(m+ncj,2,1) &
                     -(dfyy(m)+dfyy(m+ncj))*sn(m+ncj,2,2)
                tj5=(dfex(m)+dfex(m-ncj))*sn(m,2,1) &
                     +(dfey(m)+dfey(m-ncj))*sn(m,2,2) &
                     -(dfex(m)+dfex(m+ncj))*sn(m+ncj,2,1) &
                     -(dfey(m)+dfey(m+ncj))*sn(m+ncj,2,2)
!
                cndsi=sqrt(sn(m,1,1)*sn(m,1,1)+ &
                     sn(m,1,2)*sn(m,1,2)+ &
                     sn(m,1,3)*sn(m,1,3))
                cndsj=sqrt(sn(m,2,1)*sn(m,2,1)+ &
                     sn(m,2,2)*sn(m,2,2)+ &
                     sn(m,2,3)*sn(m,2,3))
                uu=(v(n,2)*sn(m,1,1)+v(n,3)*sn(m,1,2)+v(n,4)*sn(m,1,3))/ &
                     (v(n,1)*cndsi)
                vv=(v(n,2)*sn(m,2,1)+v(n,3)*sn(m,2,2)+v(n,4)*sn(m,2,3))/ &
                     (v(n,1)*cndsj)
                q2=uu**2+vv**2
                a2=cson(n)**2
                beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
                get=v(n,5)/v(n,1)             !energie totale
                ge=get-0.5*q2                 !energie interne
                coefa=coefdiag(m)-coefb(m)
                gd=(beta2-1.)/(ge*(coefa+beta2*coefb(m)))
!
                dw11=d2w1(m)+0.5*(ti1+tj1) + coefv(1,m    )*d(n-nci,1) &
                     + coefv(1,m+nci)*d(n+nci,1) &
                     + coefv(2,m    )*d(n-ncj,1) &
                     + coefv(2,m+ncj)*d(n+ncj,1)
                dw12=d2w2(m)+0.5*(ti2+tj2) + coefv(1,m    )*d(n-nci,2) &
                     + coefv(1,m+nci)*d(n+nci,2) &
                     + coefv(2,m    )*d(n-ncj,2) &
                     + coefv(2,m+ncj)*d(n+ncj,2)
                dw13=d2w3(m)+0.5*(ti3+tj3) + coefv(1,m    )*d(n-nci,3) &
                     + coefv(1,m+nci)*d(n+nci,3) &
                     + coefv(2,m    )*d(n-ncj,3) &
                     + coefv(2,m+ncj)*d(n+ncj,3)
                dw15=d2w5(m)+0.5*(ti5+tj5) + coefv(1,m    )*d(n-nci,5) &
                     + coefv(1,m+nci)*d(n+nci,5) &
                     + coefv(2,m    )*d(n-ncj,5) &
                     + coefv(2,m+ncj)*d(n+ncj,5)
                dw21=coefe(1,m    )*d(n-nci,1) &
                     +coefe(1,m+nci)*d(n+nci,1) &
                     +coefe(2,m    )*d(n-ncj,1) &
                     +coefe(2,m+ncj)*d(n+ncj,1)
                dw22=coefe(1,m    )*d(n-nci,2) &
                     +coefe(1,m+nci)*d(n+nci,2) &
                     +coefe(2,m    )*d(n-ncj,2) &
                     +coefe(2,m+ncj)*d(n+ncj,2)
                dw23=coefe(1,m    )*d(n-nci,3) &
                     +coefe(1,m+nci)*d(n+nci,3) &
                     +coefe(2,m    )*d(n-ncj,3) &
                     +coefe(2,m+ncj)*d(n+ncj,3)
                dw25=coefe(1,m    )*d(n-nci,5) &
                     +coefe(1,m+nci)*d(n+nci,5) &
                     +coefe(2,m    )*d(n-ncj,5) &
                     +coefe(2,m+ncj)*d(n+ncj,5)
!
                precon=gd*(0.5*q2*(coefa*dw11-coefb(m)*dw21) &
                     -uu*(coefa*dw12-coefb(m)*dw22) &
                     -vv*(coefa*dw13-coefb(m)*dw23) &
                     +(coefa*dw15-coefb(m)*dw25))
                d2w1(m)=dw11 + dw21 +    precon
                d2w2(m)=dw12 + dw22 + uu*precon
                d2w3(m)=dw13 + dw23 + vv*precon
                d2w5(m)=dw15 + dw25 +get*precon
             enddo
          enddo
!!!$OMP END DO
       enddo
!
!*************************************************************************
!    Calcul de l'increment implicite
!    Actualisation des variables conservatives et des flux
!    Calcul des increments de flux
!
       do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,m,ind1,ind2,wi1,wi2,wi3,wi4,wi5,ui,vi,wi,pres,fxx,fxy,fxz,fyy,fyz,fzz,fex,fey,fez)
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
!!!$OMP END DO
       enddo
!
    enddo  !fin boucle sous-iterations
!
!*************************************************************************
!      avance d'un pas de temps des variables
!*************************************************************************
!
    do k=k1,k2m1
!!!$OMP DO PRIVATE(j,n,ind1,ind2)
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
!!!$OMP END DO
    enddo
!!!$OMP END PARALLEL

    DEALLOCATE(coefe,coefv,coefdiag,coefb,d2w1,d2w2,d2w3,d2w4,d2w5)

    return
  end subroutine implimf_prcd2
end module mod_implimf_prcd2
