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
    integer          ::    njd,inc_dir(4,3),numdir
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
    double precision ::                  wi4,                 wi5,                  ww,norm,norm1
    double precision,allocatable :: coefe(:,:),   d2w1(:),   d2w2(:),   d2w3(:),   d2w4(:)
    double precision,allocatable ::    d2w5(:)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
!
    ALLOCATE(coefe(ip00,ndir))
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
    inc_dir(1,:)=(/k2m1,k2m1,k2/)
    inc_dir(2,:)=(/j2m1,j2,j2m1/)
    inc_dir(3,:)=(/i2,i2m1,i2m1/)
    inc_dir(4,:)=(/nci,ncj,nck/)

    numdir=2
    if(equat(3:4).eq.'3d') numdir=3
    fact=0.

!      calcul instationnaire avec pas de temps dual
    if(kfmg.eq.3) fact=1.5/dt1min
!
!-----initalisation--------------------------------
!
    ind1 = indc(i1m1,j1m1,k1m1)
    ind2 = indc(i2+1,j2+1,k2+1)
    do k=1,5
      do n=ind1,ind2
         d(n,k)=0.
      enddo
    enddo
    do n=ind1,ind2
       m=n-n0c
       dfxx(m)=0.
       dfxy(m)=0.
       dfxz(m)=0.
       dfex(m)=0.
       dfyy(m)=0.
       dfyz(m)=0.
       dfey(m)=0.
       dfzz(m)=0.
       dfez(m)=0.
       coefe(m,1)=0.
       coefe(m,2)=0.
       coefe(m,3)=0.
    enddo
!
!------coef diagonal ------------------------------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             coefdiag(m)=vol(n)*(fact+1./dt(n))
          enddo
       enddo
    enddo
!
!-----remplissage du coefficient diagonal par direction---------------
!
do kdir=1,numdir
   ninc=inc_dir(4,kdir)
!
    do k=k1,inc_dir(1,kdir)
       do j=j1,inc_dir(2,kdir)
          ind1 = indc(i1,j,k)
          ind2 = indc(inc_dir(3,kdir),j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             norm=1./v(n,1) ; norm1=1./v(n-ninc,1)
             uu=0.5*(v(n,2)*norm+v(n-ninc,2)*norm1)
             vv=0.5*(v(n,3)*norm+v(n-ninc,3)*norm1)
             ww=0.5*(v(n,4)*norm+v(n-ninc,4)*norm1)
             cc=0.5*(cson(n)+cson(n-ninc))
             coefe(m,kdir)=0.5*(abs(  uu*sn(m,kdir,1) &
                                    + vv*sn(m,kdir,2) &
                  + ww*sn(m,kdir,3))) + sqrt(cnds)*cc
          enddo
       enddo
    enddo
!
    do k=k1,inc_dir(1,kdir)
       do j=j1,inc_dir(2,kdir)
          ind1 = indc(i1  ,j,k)
             ind2 = indc(inc_dir(3,kdir),j,k)
          do n=ind1,ind2
             m=n-n0c
               coefdiag(m)=coefdiag(m) + coefe(m,kdir) + coefe(m+ninc,kdir)
          enddo
       enddo
    enddo
enddo
!
!*************************************************************************
!c    boucle sur les sous-iterations
!*************************************************************************
!
    do ls=1,lmx
!
!-----residu explicite------------------------------------------
!
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

       if(ityprk.ne.0) then
          do k=k1,k2m1
             do j=j1,j2m1
                ind1 = indc(i1  ,j,k)
                ind2 = indc(i2m1,j,k)
                do n=ind1,ind2
                   m=n-n0c
                   d2w1(m)=d2w1(m)-ff(n,1)
                   d2w2(m)=d2w2(m)-ff(n,2)
                   d2w3(m)=d2w3(m)-ff(n,3)
                   d2w4(m)=d2w4(m)-ff(n,4)
                   d2w5(m)=d2w5(m)-ff(n,5)
                enddo
             enddo
          enddo
       endif
!
!------direction i------------------------------------------
!
!
    do kdir=1,numdir
       do k=k1,inc_dir(1,kdir)
          do j=j1,inc_dir(2,kdir)
             ninc=inc_dir(4,kdir)
             ind1 = indc(i1,j,k)
             ind2 = indc(inc_dir(3,kdir),j,k)
             do n=ind1,ind2
                m=n-n0c
                tn1=0.5*((d(n,2)+d(n-ninc,2))*sn(m,kdir,1)    &
                     +   (d(n,3)+d(n-ninc,3))*sn(m,kdir,2)    &
                     +   (d(n,4)+d(n-ninc,4))*sn(m,kdir,3))
                tn2=0.5*((dfxx(m)+dfxx(m-ninc))*sn(m,kdir,1)  &
                     +   (dfxy(m)+dfxy(m-ninc))*sn(m,kdir,2)  &
                     +   (dfxz(m)+dfxz(m-ninc))*sn(m,kdir,3))
                tn3=0.5*((dfxy(m)+dfxy(m-ninc))*sn(m,kdir,1)  &
                     +   (dfyy(m)+dfyy(m-ninc))*sn(m,kdir,2)  &
                     +   (dfyz(m)+dfyz(m-ninc))*sn(m,kdir,3))
                tn4=0.5*((dfxz(m)+dfxz(m-ninc))*sn(m,kdir,1)  &
                     +   (dfyz(m)+dfyz(m-ninc))*sn(m,kdir,2)  &
                     +   (dfzz(m)+dfzz(m-ninc))*sn(m,kdir,3))
                tn5=0.5*((dfex(m)+dfex(m-ninc))*sn(m,kdir,1)  &
                     +   (dfey(m)+dfey(m-ninc))*sn(m,kdir,2)  &
                     +   (dfez(m)+dfez(m-ninc))*sn(m,kdir,3))
                d2w1(m)=d2w1(m) + tn1 &
                     + coefe(m,kdir)*d(n-ninc,1) &
                     + coefe(m+ninc,kdir)*d(n+ninc,1)
                d2w2(m)=d2w2(m) + tn2 &
                     + coefe(m,kdir)*d(n-ninc,2) &
                     + coefe(m+ninc,kdir)*d(n+ninc,2)
                d2w3(m)=d2w3(m) + tn3 &
                     + coefe(m,kdir)*d(n-ninc,3) &
                     + coefe(m+ninc,kdir)*d(n+ninc,3)
                d2w4(m)=d2w4(m) + tn4 &
                     + coefe(m,kdir)*d(n-ninc,4) &
                     + coefe(m+ninc,kdir)*d(n+ninc,4)
                d2w5(m)=d2w5(m) + tn5 &
                     + coefe(m,kdir)*d(n-ninc,5) &
                     + coefe(m+ninc,kdir)*d(n+ninc,5)
                d2w1(m-ninc)=d2w1(m-ninc) - tn1
                d2w2(m-ninc)=d2w2(m-ninc) - tn2
                d2w3(m-ninc)=d2w3(m-ninc) - tn3
                d2w4(m-ninc)=d2w4(m-ninc) - tn4
                d2w5(m-ninc)=d2w5(m-ninc) - tn5
             enddo
          enddo
       enddo
enddo

!
!*******************************************************************************
!
!c    calcul de l'increment implicite
!c    actualisation des variables conservatives et des flux
!c    calcul des increments de flux
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
                wi1=1./(v(n,1)+d(n,1))
                ui=v(n,2)+d(n,2)
                vi=v(n,3)+d(n,3)
                wi=v(n,4)+d(n,4)
                wi5=v(n,5)+d(n,5)

!         pression donnee par loi d'etat
                pres=gam1*(wi5-.5*wi1*(ui**2+vi**2+wi**2)-pinfl)
!
                fixx=wi1*ui**2+pres
                fixy=wi1*ui*vi
                fixz=wi1*ui*wi
                fiyy=wi1*vi**2+pres
                fiyz=wi1*vi*wi
                fizz=wi1*wi**2+pres
                fiex=wi1*ui*(wi5+pres-pinfl)
                fiey=wi1*vi*(wi5+pres-pinfl)
                fiez=wi1*wi*(wi5+pres-pinfl)
!
                norm=1./v(n,1)
                fxx=v(n,2)*(v(n,2)*norm)+ps(n)
                fxy=v(n,3)*(v(n,2)*norm)
                fxz=v(n,4)*(v(n,2)*norm)
                fyy=v(n,3)*(v(n,3)*norm)+ps(n)
                fyz=v(n,4)*(v(n,3)*norm)
                fzz=v(n,4)*(v(n,4)*norm)+ps(n)
                fex=(v(n,5)+ps(n)-pinfl)*v(n,2)*norm
                fey=(v(n,5)+ps(n)-pinfl)*v(n,3)*norm
                fez=(v(n,5)+ps(n)-pinfl)*v(n,4)*norm
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
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             v(n,:)=v(n,:)+d(n,:)
          enddo
       enddo
    enddo

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
