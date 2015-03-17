module mod_met_klnmut
  implicit none
contains
  subroutine met_klnmut( &
       l, &
       sn, &
       vol,v,mu,mut,dist, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       tprod,t, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C :  Eric Goncalves / SINUMEF
!
!     ACT
!_A   Modele de Smith (AIAA paper 97-1959)
!_A   Viscosite turbulente dynamique avec correction de non equilibre
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    dist       : arg real(ip12      ) ; distance a la paroi
!
!     LOC
!_L    t          : arg real(ip00     )  ; variable de travail
!_L    tprod      : arg real(ip00     )  ; production de k
!_L    dvxx       : arg real(ip00     )  ; grad(V)  vx,x
!_L    dvxy       : arg real(ip00     )  ; grad(V)  vx,y
!_L    dvxz       : arg real(ip00     )  ; grad(V)  vx,z
!_L    dvyx       : arg real(ip00     )  ; grad(V)  vy,x
!_L    dvyy       : arg real(ip00     )  ; grad(V)  vy,y
!_L    dvyz       : arg real(ip00     )  ; grad(V)  vy,z
!_L    dvzx       : arg real(ip00     )  ; grad(V)  vz,x
!_L    dvzy       : arg real(ip00     )  ; grad(V)  vz,y
!_L    dvzz       : arg real(ip00     )  ; grad(V)  vz,z
!
!      I/O
!_/    mut        : arg real(ip12      ) ; viscosite turbulente
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use mod_teq_gradv
    use mod_met_prod
    implicit none
    integer          ::    i,  i1,  i2,i2m1,  id
    integer          ::  inc, ind,   j,  j1,  j2
    integer          :: j2m1,  jd,   k,  k1,  k2
    integer          :: k2m1,  kd,   l,   m,   n
    integer          ::   n0, nci, ncj, nck, nid
    integer          :: nijd, njd
    double precision ::  alpha,    c1,   c14,    c2,   c22
    double precision ::  cmui1, cmui2, cmuj1, cmuj2, cmuk1
    double precision ::  cmuk2,  csk2,  dist,  dvxx,  dvxy
    double precision ::   dvxz,  dvyx,  dvyy,  dvyz,  dvzx
    double precision ::   dvzy,  dvzz,   eps,    f1,   fmu
    double precision ::     mu,   mut, mutpr,prodpr, r2sb1
    double precision ::   rack,ralpha,sigmaa,    sn,   ss2
    double precision ::      t, tprod,     v,   vol,   xl2
    double precision ::    xxi,  xxi2
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension mut(ip12),mu(ip12),dist(ip12),vol(ip11)
    dimension sn(ip31*ndir)
    dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00), &
         t(ip00),tprod(ip00)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
         cmuk1(ip21),cmuk2(ip21)
!
    ind(i,j,k) =n0 +1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    inc(id,jd,kd)=id+jd*nid+kd*nijd
!     ----------------------------------------------------------
!com  teq_gradv --> grad(v) aux points interieurs au domaine
!
!     | dvxx,dvxy,dvxz |    | du/dx du/dy du/dz |
!     | dvyx,dvyy,dvyz | => | dv/dx dv/dy dv/dz |
!     | dvzx,dvzy,dvzz |    | dw/dx dw/dy dw/dz |
!
    call    teq_gradv( &
         l, &
         sn, &
         vol,v, &
         t , &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!     ----------------------------------------------------------
!com  met_prod --> calcul de vxproduction de k
!
    call met_prod( &
         l,0, &
         mut,v, &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         tprod)
!
    n0=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
!
!     notation pour les constantes
!
!     B1 <-> cklb1
!     E2 <-> ckle2
!     sigma_k <-> sigmak
!     sigma_l <-> sigmal
!
    csk2  =2.*sqrt(2.)/cklb1
    c1=23.5
    c2=2.
    r2sb1 =sqrt(2.)/cklb1**(1./3.)
    c14=c1**4
    c22=c2**2
!
!     mutpr  : mu_t equilibre
!     prodpr : production de k avec mu_t equilibre
!     eps    : dissipation
!
    do k=k1,k2m1
       do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
             n=n+nci
             m=n-n0
             rack  =sqrt(v(n,6)/v(n,1))
             xxi   =r2sb1*rack*v(n,7)/mu(n)
             xxi2  =xxi*xxi
             xl2   =(v(n,7)/(xkappa*v(n,1)*dist(n)))**2
             f1    =exp(-50.*xl2)
             ss2   =(c22+xxi2)*xxi2
             fmu   =((c14*f1+ss2)/(c14+ss2))**0.25
             mutpr =fmu*mu(n)*xxi
             prodpr=min(tprod(m)*mutpr/mut(n),0.)
             eps   =csk2*sqrt(v(n,1)*v(n,6))*v(n,6)/v(n,7)
             alpha = prodpr/eps
             ralpha=sqrt(abs(alpha))
             sigmaa=(alpha-0.25*ralpha+0.875)/(ralpha*alpha+0.625)
             mut(n)=sigmaa*mutpr
          enddo
       enddo
    enddo
!
    return
  end subroutine met_klnmut



end module mod_met_klnmut
