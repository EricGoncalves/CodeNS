module mod_inivec
implicit none
contains
      subroutine inivec( &
                 dt,v,mu,mut, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 sn, &
                 vol, &
                 ptdual,vdual,vdual1,vdual2, &
                 cvi,cvj,cvk, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
                 pression,ztemp,cson, &
                 tnte1,tnte3,tnte4, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10)
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des tableaux pour vectorisation.
!
!     INP
!_I    reelmx     : com real             ; nombre reel grand
!_I    reelmn     : com real             ; nombre reel petit
!
!     OUT
!_O    dt         : arg real(ip11      ) ; pas de temps
!_O    v          : arg real(ip11,ip60 ) ; variables
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!_O    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_O                                        contraintes visqueuses
!_O    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_O                                        contraintes visqueuses
!_O    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_O                                        contraintes visqueuses
!_O    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_O                                        contraintes visqueuses
!_O    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_O                                        contraintes visqueuses
!_O    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_O                                        contraintes visqueuses
!_O    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_O    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_O    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!_O    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_O                                        norme egale a la surface de celle-ci
!_O    vol        : arg real(ip11      ) ; volume d'une cellule
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use constantes
      use chainecarac
implicit none
double precision :: dt
double precision :: v
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: qcx
double precision :: qcy
double precision :: qcz
double precision :: sn
double precision :: vol
double precision :: ptdual
double precision :: vdual
double precision :: vdual1
double precision :: vdual2
double precision :: cvi
double precision :: cvj
double precision :: cvk
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
double precision :: pression
double precision :: ztemp
double precision :: cson
double precision :: tnte1
double precision :: tnte3
double precision :: tnte4
double precision :: tn1
double precision :: tn2
double precision :: tn3
double precision :: tn4
double precision :: tn5
double precision :: tn6
double precision :: tn7
double precision :: tn8
double precision :: tn9
double precision :: tn10
integer :: i
integer :: m
integer :: n
!
!-----------------------------------------------------------------------
!
      real mu,mut
      dimension v(ip11,ip60)
      dimension vol(ip11),dt(ip11),pression(ip11),ztemp(ip11),cson(ip11)
      dimension mu(ip12),mut(ip12)
      dimension toxx(ip12),toxy(ip12),toxz(ip12),toyy(ip12),toyz(ip12), &
                tozz(ip12),qcx(ip12),qcy(ip12),qcz(ip12)
      dimension sn(ip31*ndir)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00), &
                tn6(ip00),tn7(ip00),tn8(ip00),tn9(ip00),tn10(ip00)
      dimension tnte1(ip11,ip60),tnte3(ip11,ip60),tnte4(ip11,ip60)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60), &
                ptdual(ip11,ip60)
      dimension cvi(ip21),cvj(ip21),cvk(ip21), &
                cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      do n=1,ip11
       v(n,1)=reelmn
       v(n,2)=reelmn
       v(n,3)=reelmn
       v(n,4)=reelmn
       v(n,5)=reelmn
       v(n,6)=reelmn
       v(n,7)=reelmn
       dt(n) =reelmx
       tnte1(n,1)=reelmx
       tnte1(n,2)=reelmx
       tnte1(n,3)=reelmx
       tnte1(n,4)=reelmx
       tnte1(n,5)=reelmx
       tnte3(n,1)=reelmx
       tnte3(n,2)=reelmx
       tnte3(n,3)=reelmx
       tnte3(n,4)=reelmx
       tnte3(n,5)=reelmx
       tnte4(n,1)=reelmx
       tnte4(n,2)=reelmx
       tnte4(n,3)=reelmx
       tnte4(n,4)=reelmx
       tnte4(n,5)=reelmx
       pression(n)=reelmn
       ztemp(n)=reelmn
       cson(n)=reelmn
      enddo
!
      do n=1,ip31*ndir
       sn(n)=reelmx
      enddo
!
      do n=1,ip11
       vol(n)=0.
      enddo
!
      do n=1,ip12
       toxx(n)=0.
       toyy(n)=0.
       tozz(n)=0.
       toxy(n)=0.
       toxz(n)=0.
       toyz(n)=0.
       qcx (n)=0.
       qcy (n)=0.
       qcz (n)=0.
       mu(n)  =0.
       mut(n) =0.
      enddo
!
      do n=1,ip21
        cvi(n)=1.
        cvj(n)=1.
        cvk(n)=1.
        cmui1(n)=1.
        cmui2(n)=1.
        cmuj1(n)=1.
        cmuj2(n)=1.
        cmuk1(n)=1.
        cmuk2(n)=1.
      enddo
!
      do n=1,ip11
       do i=1,ip60
        vdual(n,i)  = reelmx
        vdual1(n,i) = reelmx
        vdual2(n,i) = reelmx
        ptdual(n,i) = 0.
       enddo
      enddo
!
      do m= 1,ip00
       tn1 (m)=reelmx
       tn2 (m)=reelmx
       tn3 (m)=reelmx
       tn4 (m)=reelmx
       tn5 (m)=reelmx
       tn6 (m)=reelmx
       tn7 (m)=reelmx
       tn8 (m)=reelmx
       tn9 (m)=reelmx
       tn10(m)=reelmx
      enddo
!
      return
      end subroutine
end module
