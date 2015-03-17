module mod_met_iniko
  implicit none
contains
  subroutine met_iniko( &
       l,ncin,ncbd, &
       v,mut,mu,dist,mnpar, &
       sn,vol,s, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A   Initialisation de k-omega pour les modeles de WILCOX et MENTER
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use modeleturb
    use mod_met_cutked
    use mod_met_brko
    use mod_met_kocmut
    use mod_met_komut
    implicit none
    integer          ::     l,mnpar, ncbd, ncin
    double precision :: cmui1,cmui2,cmuj1,cmuj2,cmuk1
    double precision :: cmuk2, dist, dvxx, dvxy, dvxz
    double precision ::  dvyx, dvyy, dvyz, dvzx, dvzy
    double precision ::  dvzz,   mu,  mut,    s,   sn
    double precision ::     v,  vol
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension mut(ip12),mu(ip12),dist(ip12),mnpar(ip12)
    dimension ncin(ip41),ncbd(ip41),vol(ip11)
    dimension sn(ip31*ndir)
    dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00),s(ip00)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
         cmuk1(ip21),cmuk2(ip21)
!
!com  calcul de grad(V)
    call teq_gradv( &
         l, &
         sn, &
         vol,v, &
         s , &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!com  initialisation de k et omega a partir de mut
    call met_brko( &
         l, &
         mut,v, &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
!
!com  troncature des variables k-omega
    call met_cutked(l,v)
!
!com  calcul de mut a partir de k et omega
    if(kcmut.eq.4) then
       call met_komut( &
            l, &
            sn,vol,s, &
            dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
            dist,v,mu,mut, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
    else if(kcmut.eq.8) then
!         call met_kokmut(
!     &
!     &           l,
!     &           sn,vol,s,
!     &           dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz,
!     &           dist,v,mu,mut,
!     &           cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
       call met_kocmut( &
            l, &
            sn,vol,s, &
            dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
            dist,v,mu,mut, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
    endif
!
    return
  end subroutine met_iniko
end module mod_met_iniko
