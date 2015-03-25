module mod_met_inikl
  implicit none
contains
  subroutine met_inikl( &
       l,v,mut,mu,dist, &
       sn,vol,s, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A   Initialisation de k-l pour le modele de Smith
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use mod_teq_gradv
    use mod_met_cutke
    use mod_met_brkl
    use mod_met_klmut
    implicit none
  integer          :: l
  double precision ::   cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21),  cmuk1(ip21)
  double precision ::   cmuk2(ip21),   dist(ip12),   dvxx(ip00),   dvxy(ip00),   dvxz(ip00)
  double precision ::    dvyx(ip00),   dvyy(ip00),   dvyz(ip00),   dvzx(ip00),   dvzy(ip00)
  double precision ::    dvzz(ip00),     mu(ip12),    mut(ip12),      s(ip00),sn(ip31*ndir)
  double precision ::  v(ip11,ip60),    vol(ip11)
!
!-----------------------------------------------------------------------
!
    call teq_gradv( &
         l, &
         sn, &
         vol,v, &
         s, &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
    call met_brkl( &
         l, &
         mu,mut,v, &
         dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
!
!     troncature des variables k-l
!
    call met_cutke(l,v)
!
!     calcul de mut
    call met_klmut( &
         l, &
         v,mu,mut,dist)
!
    return
  end subroutine met_inikl
end module mod_met_inikl
