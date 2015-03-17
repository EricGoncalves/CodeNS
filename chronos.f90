module mod_chronos
implicit none
contains
      subroutine chronos( &
                 l,etal,mu,mut, &
                 t, &
                 dt,equat, &
                 sn,lgsnlt, &
                 vol,  &
                 sfsi,sfsj,sfsk,dism, &
                 cson)
!
!***********************************************************************
!
!_DA  DATE_C : mai 2006 -- Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul du pas de temps pour un domaine structure.
!_C    Pas de temps = min(pas de temps Euler, pas de temps visqueux).
!_C    Calcul du pas de temps local ( critere de stabilite cfl ).
!_C    Le pas de temps est calcule par point.
!
!_I    l          : arg int              ; numero de domaine
!_I    etal       : arg real             ; nombre de CFL
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    t          : arg real(ip11,ip60 ) ; variables de calcul
!_I    equat      : arg char             ; type d'equations modelisant l'ecoulement
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
      use maillage
      use proprieteflu
implicit none
integer :: inc
integer :: indc
integer :: indn
integer :: l
double precision :: etal
double precision :: t
double precision :: dt
double precision :: sn
integer :: lgsnlt
double precision :: vol
double precision :: sfsi
double precision :: sfsj
double precision :: sfsk
double precision :: dism
double precision :: cson
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
double precision :: dte
double precision :: dtv
integer :: i1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: m
integer :: mc
integer :: mn
integer :: n
integer :: n0c
integer :: n0n
integer :: n1
integer :: n2
integer :: nc
integer :: nci
integer :: ncj
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: q
double precision :: qq
double precision :: ro
double precision :: rv
double precision :: simax
double precision :: sjmax
double precision :: skmax
double precision :: smoy
!
!-----------------------------------------------------------------------
!
      double precision mu,mut,rv1
      character(len=7 ) :: equat
!
      dimension t(ip11,ip60)
      dimension mu(ip12),mut(ip12)
      dimension sn(lgsnlt,nind,ndir)
      dimension vol(ip11),cson(ip11),dt(ip11)
      dimension sfsi(ip00),sfsj(ip00),sfsk(ip00),dism(ip00)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      indn(i,j,k)=n0n+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
      inc(id,jd,kd)=id+jd*nid+kd*nijd

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
      do k=k1,k2m1
       do j=j1,j2m1
        do i=i1,i2m1
         nc=indc(i,j,k)
         mc=nc-n0c
         ro=t(nc,1)
         qq=(t(nc,2)**2+t(nc,3)**2+t(nc,4)**2)/ro**2
         q=sqrt(qq)
         dte=dism(mc)/(q+cson(nc))
         dt(nc)=etal*dte
        enddo
       enddo
      enddo

      if (equat(1:2).eq.'ns') then
       do k=k1,k2m1
        do j=j1,j2m1
         do i=i1,i2m1
          nc=indc(i,j,k)
          mc=nc-n0c
          rv1=gam*(mu(nc)/pr+mut(nc)/prt)/t(nc,1)
          rv=max(rv1,4./3.*(mu(nc)+mut(nc))/t(nc,1))
          dtv=etal*0.25*dism(mc)**2/rv
          dt(nc)=min(dt(nc),dtv)
         enddo
        enddo
       enddo
      endif

      return
      end subroutine
end module
