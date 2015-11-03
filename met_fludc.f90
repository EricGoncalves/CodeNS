module mod_met_fludc
  implicit none
contains
  subroutine met_fludc( &
       l, &
       s,mu,mut, &
       t,dtdx,dtdy,dtdz, &
       fd5x,fd5y,fd5z,fd6x,fd6y,fd6z)
!
!***********************************************************************
!
!     ACT
!_A   Calcul des densites de flux dissipatifs a partir de grad(k) et grad(e)
!_A   pour "sigma_k" et "sigma_e" constants.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    s          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    cmu        : com real             ; constante C_mu
!_I    cke1       : com real             ; constante C_epsilon-1
!_I    cke2       : com real             ; constante C_epsilon-2
!_I    alfak      : com real             ; constante sigma_k
!_I    alfak      : com real             ; constante sigma_epsilon
!
!     OUT
!
!     I/O
!_/    fd5x       : arg real(ip12     )  ; comp x grad(k) puis flux diss.
!_/    fd5y       : arg real(ip12     )  ; comp y grad(k) puis flux diss.
!_/    fd5z       : arg real(ip12     )  ; comp z grad(k) puis flux diss.
!_/    fd6x       : arg real(ip12     )  ; comp x grad(e) puis flux diss.
!_/    fd6y       : arg real(ip12     )  ; comp y grad(e) puis flux diss.
!_/    fd6z       : arg real(ip12     )  ; comp z grad(e) puis flux diss.
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use chainecarac
    implicit none
    integer          ::    i,  i1,i1m1,  i2,i2m1
    integer          :: imax,imin,ind1,ind2,   j
    integer          ::   j1,j1m1,  j2,j2m1,jmax
    integer          :: jmin,   k,  k1,k1m1,  k2
    integer          :: k2m1,kmax,kmin,   l,   n
    integer          ::   n0, nid,nijd, njd
    double precision :: dtdx,dtdy,dtdz,fd5x,fd5y
    double precision :: fd5z,fd6x,fd6y,fd6z,  mu
    double precision ::  mut,   s,   t
!
!-----------------------------------------------------------------------
!
    dimension s(ip11,ip60)
    dimension mut(ip12),mu(ip12), &
         fd5x(ip12),fd5y(ip12),fd5z(ip12), &
         fd6x(ip12),fd6y(ip12),fd6z(ip12)
    dimension t(ip00),dtdx(ip00),dtdy(ip00),dtdz(ip00)
!
    n0=npc(l)
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    i1m1=i1-1
    j1m1=j1-1
    k1m1=k1-1
!
    i2m1=i2-1
    j2m1=j2-1
    k2m1=k2-1
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd= nid*njd
!
    imin=i1m1
    imax=i2
    jmin=j1m1
    jmax=j2
    kmin=k1m1
    kmax=k2
!
    if (equat(3:5).eq.'2di') then
       imin=i1
       imax=i2m1
    endif
    if (equat(3:5).eq.'2dj') then
       jmin=j1
       jmax=j2m1
    endif
    if (equat(3:5).eq.'2dk') then
       kmin=k1
       kmax=k2m1
    endif
!
    ind1=indc(i1  ,j1  ,k1  )
    ind2=indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       fd5x(n)=(mu(n)+mut(n)/alfak)*fd5x(n)
       fd5y(n)=(mu(n)+mut(n)/alfak)*fd5y(n)
       fd5z(n)=(mu(n)+mut(n)/alfak)*fd5z(n)
    enddo
!
    ind1=indc(i1  ,j1  ,k1  )
    ind2=indc(i2m1,j2m1,k2m1)
    do n=ind1,ind2
       fd6x(n)=(mu(n)+mut(n)/alfae)*fd6x(n)
       fd6y(n)=(mu(n)+mut(n)/alfae)*fd6y(n)
       fd6z(n)=(mu(n)+mut(n)/alfae)*fd6z(n)
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine met_fludc
end module mod_met_fludc
