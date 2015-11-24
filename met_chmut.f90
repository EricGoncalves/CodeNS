module mod_met_chmut
  implicit none
contains
  subroutine met_chmut( &
       l, &
       v,mu,mut,dist,mnpar,utau)
!
!***********************************************************************
!
!     ACT
!_A   Modele de Chien
!_A   Calcul de la viscosite turbulente a partir de k et epsilon
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use chainecarac
    implicit none
    integer          ::           i,         i1,         i2,       i2m1
    integer          ::           j,         j1,         j2,       j2m1
    integer          ::           k,         k1,         k2,       k2m1
    integer          ::           l,mnpar(ip12),         mp,          n,         n0
    integer          ::         nci,        ncj,        nck,        nid,       nijd
    integer          ::         njd
    double precision ::   dist(ip12),         fmu,    mu(ip12),   mut(ip12),  utau(ip42)
    double precision :: v(ip11,ip60),       yplus
!
!-----------------------------------------------------------------------
!
!


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
    do k=k1,k2m1
       do j=j1,j2m1
          n=ind(i1-1,j,k)
          do i=i1,i2m1
             n=n+nci
             mp    =mnpar(n)
!         yplus =utau(mp)*dist(n)*v(n,1)/mu(n)
             yplus=max(abs(utau(mp)),utaumin)*dist(n)*v(n,1)/mu(n)
             fmu   =1.-exp(-0.0115*yplus)
             mut(n)=cmu*fmu*(v(n,6)**2)/v(n,7)
          enddo
       enddo
    enddo
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
    end function ind
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine met_chmut
end module mod_met_chmut
