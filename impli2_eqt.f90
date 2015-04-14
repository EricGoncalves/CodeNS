module mod_impli2_eqt
  implicit none
contains
  subroutine impli2_eqt( &
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
    integer          ::     ibalai,        id,      ind1,      ind2,         j
    integer          ::         j1,      j1m1,        j2,      j2m1,        jd
    integer          ::          k,        k1,      k1m1,        k2,      k2m1
    integer          ::         kd,      kdir,         l,      ldom,    lgsnlt
    integer          ::         li,        lj,         m,        mb,        mf
    integer          ::        mfb,        mt,         n,       n0c,       nci
    integer          :: ncin(ip41),       ncj,       nck,      ncyc,        ni
    integer          ::        nid,      nijd,      ninc,       njd,        no
    double precision ::                   ai,         alpha(ip00),         beta6(ip00),         beta7(ip00),                  bi
    double precision ::                  cci,          cfke(ip13),                 cmt,                cnds,      coefdiag(ip00)
    double precision ::                  di6,                 di7,                 dj6,                 dj7,            dt(ip11)
    double precision ::           dwi6(ip00),          dwi7(ip00),                fact,            mu(ip12),           mut(ip12)
    double precision ::             rv(ip00),sn(lgsnlt,nind,ndir),                  td,                 tmi,                 tmj
    double precision ::                  tpi,                 tpj,        u(ip11,ip60),            u1(ip00),            u2(ip00)
    double precision ::             u3(ip00),                  uu,        v(ip11,ip60),           vol(ip11),                  vv
    double precision ::                   ww
    double precision,allocatable :: coefe(:,:)
!
!-----------------------------------------------------------------
!
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

!     nombre de balayage par direction
    ibalai=2

!     constante instationnaire dts
    fact=1.5
!     fact=11./6.  !ordre 3

!     constante du modele
    if(equatt(1:3).eq.'2JL') then
       cmt=1.
    elseif(equatt(1:3).eq.'2Sm') then
       cmt=1.43
    elseif(equatt(1:3).eq.'2WL') then
       cmt=2.
    elseif(equatt(1:3).eq.'2MT') then
       cmt=1.17
    elseif(equatt(1:3).eq.'1SA') then
       cmt=1.
    elseif(equatt(1:3).eq.'2KO') then
       cmt=2.
    endif
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
    do k=k1,k2m1
       do j=j1,j2m1
          ind1 = indc(i1  ,j,k)
          ind2 = indc(i2m1,j,k)
          do n=ind1,ind2
             m=n-n0c
             coefdiag(m)=coefdiag(m) + coefe(1,m) + coefe(1,m+nci) &
                  + coefe(2,m) + coefe(2,m+ncj)
          enddo
       enddo
    enddo
!
!------calcul instationnaire avec dts-----------------------------
!
    if(kfmg.eq.3) then
       do k=k1,k2m1
          do j=j1,j2m1
             ind1 = indc(i1  ,j,k)
             ind2 = indc(i2m1,j,k)
             do n=ind1,ind2
                m=n-n0c
                coefdiag(m)=coefdiag(m) + fact*vol(n)/dt1min
             enddo
          enddo
       enddo
    endif
!
!*******************************************************************
!                          CAS 2D
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
                td=-u2(m+ncj)+u2(m)
                tpj=coefe(2,m+ncj)-u2(m+ncj)
                tmj=coefe(2,m)+u2(m)
                di6=-u(n,6)+td*dwi6(m)+tmj*dwi6(m-ncj)+tpj*dwi6(m+ncj)
                di7=-u(n,7)+td*dwi7(m)+tmj*dwi7(m-ncj)+tpj*dwi7(m+ncj)
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
                td=-u1(m+nci)+u1(m)
                tpi=coefe(1,m+nci)-u1(m+nci)
                tmi=coefe(1,m)+u1(m)
                dj6=-u(n,6)+td*dwi6(m)+tmi*dwi6(m-nci)+tpi*dwi6(m+nci)
                dj7=-u(n,7)+td*dwi7(m)+tmi*dwi7(m-nci)+tpi*dwi7(m+nci)
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
!-----lois de paroi------------------------------------------------
!
    if(lparoi.eq.1) then
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
  end subroutine impli2_eqt
end module mod_impli2_eqt
