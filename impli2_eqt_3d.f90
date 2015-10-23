module mod_impli2_eqt_3d
  implicit none
contains
  subroutine impli2_eqt_3d( &
       l,u,dt,v, &
       mu,mut,cfke,ncin, &
       ncyc, &
       sn,lgsnlt, &
       vol, &
       dwi6,dwi7,u1,u2,u3, &
       rv,coefdiag,alpha,beta6,beta7)
!
!******************************************************************
!
!_DATE  juin 2004 - Eric GONCALVES / LEGI
!
!     Phase implicite sans matrice avec relaxation de type
!     Jacobi par lignes alternees.
!     Version 3d pour calculs paralelles
!
!-----parameters figes-------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use schemanum
    use chainecarac
    use modeleturb
    implicit none
    integer          ::          i,        i1,      i1m1,        i2,      i2m1
    integer          ::     ibalai,      ind1,      ind2,         j
    integer          ::         j1,      j1m1,        j2,      j2m1
    integer          ::          k,        k1,      k1m1,        k2,      k2m1
    integer          ::       kdir,         l,      ldom,    lgsnlt
    integer          ::         li,        lj,        lk,         m,        mb
    integer          ::         mf,       mfb,        mt,         n,       n0c
    integer          ::        nci,ncin(ip41),       ncj,       nck,      ncyc
    integer          ::         ni,       nid,      nijd,      ninc,       njd
    integer          ::         no
    double precision ::                   ai,         alpha(ip00),         beta6(ip00),         beta7(ip00),                  bi
    double precision ::                  cci,          cfke(ip13),                 cmt,                cnds,      coefdiag(ip00)
    double precision ::                  di6,                 di7,                 dj6,                 dj7,                 dk6
    double precision ::                  dk7,            dt(ip11),          dwi6(ip00),          dwi7(ip00),                fact
    double precision ::             mu(ip12),           mut(ip12),            rv(ip00),sn(lgsnlt,nind,ndir),                  td
    double precision ::                  tmi,                 tmj,                 tmk,                 tpi,                 tpj
    double precision ::                  tpk,        u(ip11,ip60),            u1(ip00),            u2(ip00),            u3(ip00)
    double precision ::                   uu,        v(ip11,ip60),           vol(ip11),                  vv,                  ww
    double precision,allocatable :: coefe(:,:)
!
!-----------------------------------------------------------------
!
!
!



    ALLOCATE(coefe(ndir,ip00))

    n0c=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    nid=id2(l)-id1(l)+1
    njd=jd2(l)-jd1(l)+1
    nijd=nid*njd
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    nci=inc(1,0,0)
    ncj=inc(0,1,0)
    nck=inc(0,0,1)
!
!     nombre de balayage par direction
    ibalai=2
!
!     constante du modele
!
    select case(equatt(1:3))
    case('2JL')
       cmt=1.
    case('2Sm')
       cmt=1.43
    case('2WL')
       cmt=2.
    case('2MT')
       cmt=1.17
    case('1SA')
       cmt=1.
    case('2KO')
       cmt=2.
    end select
!
!-----initialisation-----------------------------------------------
!
    ind1=indc(i1m1,j1m1,k1m1)
    ind2=indc(i2+1,j2+1,k2+1)
    do n=ind1,ind2
       m=n-n0c
       dwi6(m)=0.
       dwi7(m)=0.
       coefe(1,m)=0.
       coefe(2,m)=0.
       coefe(3,m)=0.
       rv(m)=0.
       u1(m)=0.
       u2(m)=0.
       u3(m)=0.
       alpha(m)=0.
       beta6(m)=0.
       beta7(m)=0.
    enddo
!
!----calculs du rayon spectral visqueux-------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             rv(m)=(mu(n)+mut(n)/cmt)/v(n,1)
!         cfke(n)=min(cfke(n),100.)
             coefdiag(m)=vol(n)/dt(n) + cfke(n)*vol(n)
          enddo
       enddo
    enddo
!
!*****************************************************************
!-----remplissage des coefficients par direction
!*****************************************************************
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
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu=0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv=0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww=0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             u1(m)=0.5*(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))
             coefe(kdir,m)=abs(u1(m)) &
                  + (rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
          enddo
       enddo
    enddo
!
    kdir=2
    ninc=ncj
!
    do k=k1,k2m1
       do j=j1,j2
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             cnds=sn(m,kdir,1)*sn(m,kdir,1)+ &
                  sn(m,kdir,2)*sn(m,kdir,2)+ &
                  sn(m,kdir,3)*sn(m,kdir,3)
             uu  = 0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv  = 0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww  = 0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             u2(m)=0.5*(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))
             coefe(kdir,m)=abs(u2(m)) &
                  + (rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
          enddo
       enddo
    enddo
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
             uu  = 0.5*(v(n,2)/v(n,1)+v(n-ninc,2)/v(n-ninc,1))
             vv  = 0.5*(v(n,3)/v(n,1)+v(n-ninc,3)/v(n-ninc,1))
             ww  = 0.5*(v(n,4)/v(n,1)+v(n-ninc,4)/v(n-ninc,1))
             u3(m)=0.5*(uu*sn(m,kdir,1)+vv*sn(m,kdir,2)+ww*sn(m,kdir,3))
             coefe(kdir,m)=abs(u3(m)) &
                  + (rv(m)+rv(m-ninc))*cnds/(vol(n)+vol(n-ninc))
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
             coefdiag(m)=coefdiag(m) + coefe(1,m) + coefe(1,m+nci) &
                  + coefe(2,m) + coefe(2,m+ncj) &
                  + coefe(3,m) + coefe(3,m+nck)
          enddo
       enddo
    enddo
!
!c-----calcul instationnaire avec dts-----------------------------
!
!       if(kfmg.eq.3) then
!        fact=1.5
!        fact=11./6.  !ordre 3
!        do k=k1,k2m1
!         do j=j1,j2m1
!          ind1 = indc(i1  ,j,k)
!          ind2 = indc(i2m1,j,k)
!          do n=ind1,ind2
!           m=n-n0c
!           coefdiag(m)=coefdiag(m) + fact*vol(n)/dt1min
!          enddo
!         enddo
!        enddo
!       endif
!
!*******************************************************************
!                          CAS 3D
!     inversion du systeme par direction  - algorithme de Thomas
!*******************************************************************
!
!-----inversion direction i-----------------------------------------
!
    do li=1,ibalai
!
       do k=k1,k2m1
          do i=i2m1,i1,-1
             ind1=indc(i,j1  ,k)
             ind2=indc(i,j2m1,k)
             do n=ind1,ind2,ncj
                m=n-n0c
                td=-u2(m+ncj)+u2(m)-u3(m+nck)+u3(m)
                tpj=coefe(2,m+ncj)-u2(m+ncj)
                tmj=coefe(2,m)+u2(m)
                tpk=coefe(3,m+nck)-u3(m+nck)
                tmk=coefe(3,m)+u3(m)
                di6=-u(n,6)+td*dwi6(m)+tmj*dwi6(m-ncj)+tpj*dwi6(m+ncj) &
                     +tmk*dwi6(m-nck)+tpk*dwi6(m+nck)
                di7=-u(n,7)+td*dwi7(m)+tmj*dwi7(m-ncj)+tpj*dwi7(m+ncj) &
                     +tmk*dwi7(m-nck)+tpk*dwi7(m+nck)
                ai=coefe(1,m)+u1(m)
                bi=coefdiag(m)+u1(m+nci)-u1(m)
                cci=-coefe(1,m+nci)+u1(m+nci)
                alpha(m)=ai/(bi+cci*alpha(m+nci))
                beta6(m)=(di6-cci*beta6(m+nci))/(bi+cci*alpha(m+nci))
                beta7(m)=(di7-cci*beta7(m+nci))/(bi+cci*alpha(m+nci))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do i=i1,i2m1
             ind1=indc(i,j1  ,k)
             ind2=indc(i,j2m1,k)
             do n=ind1,ind2,ncj
                m=n-n0c
                dwi6(m)=alpha(m)*dwi6(m-nci)+beta6(m)
                dwi7(m)=alpha(m)*dwi7(m-nci)+beta7(m)
             enddo
          enddo
       enddo
!
    enddo
!
!-----inversion direction j------------------------------------
!
    ind1=indc(i1m1,j1m1,k1m1)-n0c
    ind2=indc(i2+1,j2+1,k2+1)-n0c
    do m=ind1,ind2
       alpha(m)=0.
       beta6(m)=0.
       beta7(m)=0.
    enddo
!
    do lj=1,ibalai
!
       do k=k1,k2m1
          do j=j2m1,j1,-1
             ind1=indc(i1  ,j,k)
             ind2=indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                td=-u1(m+nci)+u1(m)-u3(m+nck)+u3(m)
                tpi=coefe(1,m+nci)-u1(m+nci)
                tmi=coefe(1,m)+u1(m)
                tpk=coefe(3,m+nck)-u3(m+nck)
                tmk=coefe(3,m)+u3(m)
                dj6=-u(n,6)+td*dwi6(m)+tmi*dwi6(m-nci)+tpi*dwi6(m+nci) &
                     +tmk*dwi6(m-nck)+tpk*dwi6(m+nck)
                dj7=-u(n,7)+td*dwi7(m)+tmi*dwi7(m-nci)+tpi*dwi7(m+nci) &
                     +tmk*dwi7(m-nck)+tpk*dwi7(m+nck)
                ai=coefe(2,m)+u2(m)
                bi=coefdiag(m)+u2(m+ncj)-u2(m)
                cci=-coefe(2,m+ncj)+u2(m+ncj)
                alpha(m)=ai/(bi+cci*alpha(m+ncj))
                beta6(m)=(dj6-cci*beta6(m+ncj))/(bi+cci*alpha(m+ncj))
                beta7(m)=(dj7-cci*beta7(m+ncj))/(bi+cci*alpha(m+ncj))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1  ,j,k)
             ind2=indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                dwi6(m)=alpha(m)*dwi6(m-ncj)+beta6(m)
                dwi7(m)=alpha(m)*dwi7(m-ncj)+beta7(m)
             enddo
          enddo
       enddo
!
    enddo
!
!------inversion direction k--------------------------------------
!
    ind1=indc(i1m1,j1m1,k1m1)-n0c
    ind2=indc(i2+1,j2+1,k2+1)-n0c
    do m=ind1,ind2
       alpha(m)=0.
       beta6(m)=0.
       beta7(m)=0.
    enddo
!
    do lk=1,ibalai
!
       do k=k2m1,k1,-1
          do j=j1,j2m1
             ind1=indc(i1  ,j,k)
             ind2=indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                td=-u1(m+nci)+u1(m)-u2(m+ncj)+u2(m)
                tpi=coefe(1,m+nci)-u1(m+nci)
                tmi=coefe(1,m)+u1(m)
                tpj=coefe(2,m+ncj)-u2(m+ncj)
                tmj=coefe(2,m)+u2(m)
                dk6=-u(n,6)+td*dwi6(m)+tmi*dwi6(m-nci)+tpi*dwi6(m+nci) &
                     +tmj*dwi6(m-ncj)+tpj*dwi6(m+ncj)
                dk7=-u(n,7)+td*dwi7(m)+tmi*dwi7(m-nci)+tpi*dwi7(m+nci) &
                     +tmj*dwi7(m-ncj)+tpj*dwi7(m+ncj)
                ai=coefe(3,m)+u3(m)
                bi=coefdiag(m)+u3(m+nck)-u3(m)
                cci=-coefe(3,m+nck)+u3(m+nck)
                alpha(m)=ai/(bi+cci*alpha(m+nck))
                beta6(m)=(dk6-cci*beta6(m+nck))/(bi+cci*alpha(m+nck))
                beta7(m)=(dk7-cci*beta7(m+nck))/(bi+cci*alpha(m+nck))
             enddo
          enddo
       enddo
!
       do k=k1,k2m1
          do j=j1,j2m1
             ind1=indc(i1  ,j,k)
             ind2=indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                dwi6(m)=alpha(m)*dwi6(m-nck)+beta6(m)
                dwi7(m)=alpha(m)*dwi7(m-nck)+beta7(m)
             enddo
          enddo
       enddo
!
    enddo
!
!-----lois de paroi------------------------------------------------
!
    if((lparoi.eq.1).or.(lparoi.eq.2)) then
       nbd=0
       do no=1,mtbx
          mfb=nba(no)
          ldom=ndlb(mfb)
          if((cl(mfb)(1:2).eq.'lp').and.(l.eq.ldom)) then
             nbd=nbd+1
             lbd(nbd)=mfb
          endif
       enddo
!
       do mf=1,nbd
          mfb=lbd(mf)
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             ni=ncin(mb)-n0c
             if(equatt(1:3).eq.'1SA') then
                dwi6(ni)=0.
             endif
             dwi7(ni)=0.
          enddo
       enddo
!
    endif
!
!-----avance en temps------------------------------------------------
!
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             v(n,6)=v(n,6)+dwi6(m)
             v(n,7)=v(n,7)+dwi7(m)
          enddo
       enddo
    enddo

    DEALLOCATE(coefe)

    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine impli2_eqt_3d
end module mod_impli2_eqt_3d
