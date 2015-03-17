module mod_met_ini
  implicit none
contains
  subroutine met_ini( &
       l,v,mut,mu, &
       sn,vol,s, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT      initialisation k-epsilon
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use mod_teq_gradv
    use mod_met_brad
    implicit none
    integer          :: l
    double precision :: cmui1,cmui2,cmuj1,cmuj2,cmuk1
    double precision :: cmuk2, dvxx, dvxy, dvxz, dvyx
    double precision ::  dvyy, dvyz, dvzx, dvzy, dvzz
    double precision ::    mu,  mut,    s,   sn,    v
    double precision ::   vol
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension mut(ip12),mu(ip12)
    dimension sn(ip31*ndir), &
         vol(ip11)
    dimension dvxx(ip00),dvxy(ip00),dvxz(ip00), &
         dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00), &
         s(ip00)
    dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
         cmuk1(ip21),cmuk2(ip21)
!
    call teq_gradv( &
         l, &
         sn, &
         vol,v, &
         s , &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    call met_brad( &
         l, &
         mut,v, &
         dvxy,dvxz,dvyx,dvyz,dvzx,dvzy)
!
!     troncature des variables k-epsilon
!
!      call met_cut(l,v)
!
    return
  end subroutine met_ini
end module mod_met_ini
