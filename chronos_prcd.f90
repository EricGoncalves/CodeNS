module mod_chronos_prcd
  implicit none
contains
  subroutine chronos_prcd( &
       l,etal,mu,mut, &
       t, &
       dt,equat, &
       sn,lgsnlt, &
       vol,  &
       sfsi,sfsj,sfsk,dism, &
       cson,ps)
!
!***********************************************************************
!
!_DA  DATE_C : mai 2006 -- Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul du pas de temps pour un domaine structure.
!_C    Pas de temps = min(pas de temps euler, pas de temps visqueux).
!_C    Calcul du pas de temps local (critere de stabilite cfl).
!_C    Le pas de temps est calcule par point.
!_C    Preconditionnement basse vitesse de Turkel.
!
!_I    l          : arg int              ; numero de domaine
!_I    etal       : arg real             ; nombre de CFL
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    t          : arg real(ip11,ip60 ) ; variables de calcul
!_I    equat      : arg char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    sn         : arg real(lgsnlt,
!_I                          nind,ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    lgsnlt     : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
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
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    gam1       : com real             ; rap chal spec -1
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
!
!     OUT
!_O    dt         : arg real(ip11      ) ; pas de temps
!
!     LOC
!_L    sfsi       : arg real(ip00      ) ; surface au carre d'une facette i
!_L    sfsj       : arg real(ip00      ) ; surface au carre d'une facette j
!_L    sfsk       : arg real(ip00      ) ; surface au carre d'une facette k
!_L    dism       : arg real(ip00      ) ; distance caracteristique d'une maille
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use schemanum
    use definition
    use maillage
    use proprieteflu
    implicit none
    integer          ::       i,     i1,     i2,   i2m1,     id
    integer          ::     inc,   indc,   indn,isortie,  ivisq
    integer          ::       j,     j1,     j2,   j2m1,     jd
    integer          ::       k,     k1,     k2,   k2m1,     kd
    integer          ::       l, lgsnlt,      m,     mc,     mn
    integer          ::       n,    n0c,    n0n,     n1,     n2
    integer          ::      nc,    nci,    ncj,    nck,    nid
    integer          ::    nijd,    njd
    double precision ::   Lref,    a2, beta2,beta2e, betau
    double precision ::   betv,  cflc,   cpi,  cson,  dism
    double precision ::    dpi,   dpj,   dpm,    dt,   dte
    double precision ::    dtv,  etal,    mu,   mut,    ps
    double precision ::      q,    q2,  qinf,  rlam,    rv
    double precision ::    rv1,  sfsi,  sfsj,  sfsk, simax
    double precision ::  sjmax, skmax,  smoy,    sn,     t
    double precision ::     uu,   vol,    vv,    ww,   xm0
    double precision ::  xmach,xmach2
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
    dimension t(ip11,ip60)
    dimension dt(ip11),vol(ip11),cson(ip11),ps(ip11)
    dimension mu(ip12),mut(ip12)
    dimension sn(lgsnlt,nind,ndir)
    dimension sfsi(ip00),sfsj(ip00),sfsk(ip00),dism(ip00)
!
    indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    indn(i,j,k)=n0n+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    inc(id,jd,kd)=id+jd*nid+kd*nijd
!
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     :: beta2v
    ALLOCATE(beta2v(ip21))

    n0c =npc(l)
    n0n =npn(l)
    i1  =ii1(l)
    i2  =ii2(l)
    j1  =jj1(l)
    j2  =jj2(l)
    k1  =kk1(l)
    k2  =kk2(l)
!
    nid  = id2(l)-id1(l)+1
    njd  = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    n1=indc(i1,j1,k1)
    n2=indc(i2,j2,k2)
    do n=n1,n2
       dt(n)=0.
    enddo
!
    do k=k1,k2m1
       do j=j1,j2m1
!!$OMP SIMD
          do i=i1,i2
             n=indn(i,j,k)
             m=n-n0n
             sfsi(m)=sn(m,1,1)*sn(m,1,1)+ &
                  sn(m,1,2)*sn(m,1,2)+ &
                  sn(m,1,3)*sn(m,1,3)
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       do j=j1,j2
!!$OMP SIMD
          do i=i1,i2m1
             n=indn(i,j,k)
             m=n-n0n
             sfsj(m)=sn(m,2,1)*sn(m,2,1)+ &
                  sn(m,2,2)*sn(m,2,2)+ &
                  sn(m,2,3)*sn(m,2,3)
          enddo
       enddo
    enddo
!
    do k=k1,k2
       do j=j1,j2m1
!!$OMP SIMD
          do i=i1,i2m1
             n=indn(i,j,k)
             m=n-n0n
             sfsk(m)=sn(m,3,1)*sn(m,3,1)+ &
                  sn(m,3,2)*sn(m,3,2)+ &
                  sn(m,3,3)*sn(m,3,3)
          enddo
       enddo
    enddo
!
    do k=k1,k2m1
       do j=j1,j2m1
          mn=indn(i1-1,j,k)-n0n
          mc=indc(i1-1,j,k)-n0c
!!$OMP SIMD
          do i=i1,i2m1
             mn=mn+nci
             mc=mc+nci
             nc=mc+n0c
             simax=max(sfsi(mn+nci),sfsi(mn))
             sjmax=max(sfsj(mn+ncj),sfsj(mn))
             skmax=max(sfsk(mn+nck),sfsk(mn))
             smoy =sqrt(simax+sjmax+skmax)
             dism(mc)=vol(nc)/smoy
          enddo
       enddo
    enddo
!
!-----precond avec formulation max grad P pour toutes les faces
!      do k=k1,k2m1
!       do j=j1,j2m1
!        do i=i1,i2m1
!         nc=indc(i,j,k)
!         mc=nc-n0c
!         dpi=max(abs(ps(nc+nci)-ps(nc)),abs(ps(nc-nci)-ps(nc)))
!         dpj=max(abs(ps(nc+ncj)-ps(nc)),abs(ps(nc-ncj)-ps(nc)))
!         dpm=max(dpi,dpj)
!         q2=(t(nc,2)**2+t(nc,3)**2+t(nc,4)**2)/t(nc,1)**2
!         q=sqrt(q2)
!         a2=cson(nc)**2
!         beta2=min(max(q2/a2,cte*dpm/(t(nc,1)*a2)),1.)
!         beta2v(mc)=beta2
!        enddo
!       enddo
!      enddo
!
!
    qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
    do k=k1,k2m1
       do j=j1,j2m1
          do i=i1,i2m1
             nc=indc(i,j,k)
             mc=nc-n0c
             q2=(t(nc,2)**2+t(nc,3)**2+t(nc,4)**2)/t(nc,1)**2
             q=sqrt(q2)
             a2=cson(nc)**2
             beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
!         beta2e=max(q2/a2,cte*qinf**2/a2)
!         beta2=min(max(beta2e,beta2v(m)),1.)
             rlam=0.5*((1.+beta2)*q+sqrt(((1.-beta2)*q)**2+4.*beta2*a2))
             dte=dism(mc)/rlam
             dt(nc)=etal*dte
          enddo
       enddo
    enddo
!
    if (equat(1:2).eq.'ns') then
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                nc=indc(i,j,k)
                mc=nc-n0c
                rv1=gam*(mu(nc)/pr+mut(nc)/prt)/t(nc,1)
                rv=max(rv1,4./3.*(mu(nc)+mut(nc))/t(nc,1))
                dtv=etal*0.5*dism(mc)**2/rv
                dt(nc)=min(dt(nc),dtv)
             enddo
          enddo
       enddo
    endif
!
!--------beta2 pour instationnaire---------------------------------
!
    if((kfmg.eq.3).and.(kvisq.eq.1)) then
       Lref=1.25
       cpi=3.14159265
       xm0=0.6
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                nc=indc(i,j,k)
                mc=nc-n0c
                q2=(t(nc,2)**2+t(nc,3)**2+t(nc,4)**2)/t(nc,1)**2
                q=sqrt(q2)
                xmach=q/cson(nc)
                xmach2=xmach**2
!           cflc=cson(nc)*dt(nc)/dism(mc)
!           betau=(Lref/(cpi*dt(nc)*cson(nc)))**2
!           betau=xmach/cpi**2
!           betau=xmach**2+1./cflc**2
                betau=xmach2*(1.+xmach2*(1.-xm0**2)/xm0**4)
                beta2v(mc)=min(max(beta2v(mc),betau),1.)
             enddo
          enddo
       enddo
    endif
!
!--------beta2 pour ecoulements visqueux----------------------------
!
    ivisq=0
    if(ivisq.eq.1) then
       do k=k1,k2m1
          do j=j1,j2m1
             mn=indn(i1-1,j,k)-n0n
             mc=indc(i1-1,j,k)-n0c
             do i=i1,i2m1
                mn=mn+nci
                mc=mc+nci
                nc=mc+n0c
                uu=t(nc,2)/t(nc,1)
                vv=t(nc,3)/t(nc,1)
                ww=t(nc,4)/t(nc,1)
                q2=uu**2+vv**2+ww**2
                a2=cson(nc)**2
                q=sqrt(q2)
                betv=7.
                beta2v(mc)=betv
             enddo
          enddo
       enddo
    endif
!
    isortie=0
    if(isortie.eq.1) then
       write(6,'("===>chronos_prcd")')
       k=1
!       i=175
       do i=i1,i2m1
          do j=j1,j2m1
             n=indc(i,j,k)
             m=n-n0c
             write(6,'(2i4,i6,(1pe12.4))') &
                  i,j,n,beta2v(m)
          enddo
       enddo
    endif

    DEALLOCATE(beta2v)

    return
  end subroutine chronos_prcd
end module mod_chronos_prcd
