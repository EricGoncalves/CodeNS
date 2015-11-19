module mod_atctranske
  implicit none
contains
  subroutine atctranske( &
       l,v,mu,mut,mnpar,fgam)
!
!***********************************************************************
!
!     ACT
!_A    Impose une transition fixee definie par la fonction "fgam"
!_A    en chaque point de la paroi.
!
!     VAL
!_V    Uniquement pour modeles algebriques
!
!     INP
!_I    id1        : arg int              ; indice min en i fictif
!_I    ii1        : arg int              ; indice min en i reel
!_I    ii2        : arg int              ; indice max en i reel
!_I    id2        : arg int              ; indice max en i fictif
!_I    jd1        : arg int              ; indice min en j fictif
!_I    jj1        : arg int              ; indice min en j reel
!_I    jj2        : arg int              ; indice max en j reel
!_I    jd2        : arg int              ; indice max en j fictif
!_I    kd1        : arg int              ; indice min en k fictif
!_I    kk1        : arg int              ; indice min en k reel
!_I    kk2        : arg int              ; indice max en k reel
!_I    kd2        : arg int              ; indice max en k fictif
!_I    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_I                                        stockees du point de rattach normale
!_I    fgam       : arg real(ip42      ) ; fonction d'intermittence pour transition
!
!     OUT
!_/    mut        : arg real(ip12      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    implicit none
    integer          ::          i1,         i2,       i2m1,       ind1
    integer          ::        ind2,          j,         j1,         j2,       j2m1
    integer          ::           k,         k1,         k2,       k2m1,          l
    integer          :: mnpar(ip12),       mpar,          n,        n0c,        nid
    integer          ::        nijd,        njd
    double precision ::   fgam(ip42),    mu(ip12),   mut(ip12),v(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!
    n0c=npc(l)
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
    do k=k1,k2m1
       do j=j1,j2m1
          ind1=indc(i1,j,k)
          ind2=indc(i2m1,j,k)
          do n=ind1,ind2
             mpar=mnpar(n)
             if(mpar.ge.1) then
                if(fgam(mpar).lt.1.e-3) then
                   mut(n)=1.e-3*mu(n)
                   v(n,6)=epsk
                   v(n,7)=epse
                end if
!         mut(n)=max(fgam(mpar)*mut(n),1.e-3*mu(n))
             end if
          enddo
       enddo
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function indc
  end subroutine atctranske
end module mod_atctranske
