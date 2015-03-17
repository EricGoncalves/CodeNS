module mod_cccca
  implicit none
contains
  subroutine cccca( &
       l,equat, &
       v,mut, &
       vv1,vv2,vv3,vv4,vv5,vv6,vv7,mmut)
!
!***********************************************************************
!
!     ACT
!_A    Transfert de variables aux centres a variables aux centres.
!_A    y compris pour les mailles fictives.
!
!_I    l          : arg int              ; numero de domaine
!_I    equat      : arg char             ; type d'equations modelisant l'ecoulement
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    nnc        : com int (lt        ) ; nombre de cellules du dom (dont fic.)
!
!     OUT
!_O    vv1        : arg real(ip00      ) ; variable de calcul numero 1
!_O    vv2        : arg real(ip00      ) ; variable de calcul numero 2
!_O    vv3        : arg real(ip00      ) ; variable de calcul numero 3
!_O    vv4        : arg real(ip00      ) ; variable de calcul numero 4
!_O    vv5        : arg real(ip00      ) ; variable de calcul numero 5
!_O    mmut       : arg real(ip00      ) ; viscosite turbulente
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    implicit none
    integer          ::  l, m, n,n0
    double precision :: mmut, mut,   v, vv1, vv2
    double precision ::  vv3, vv4, vv5, vv6, vv7
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: equat
    dimension v(ip11,ip60),mut(ip12)
    dimension vv1(ip00),vv2(ip00),vv3(ip00),vv4(ip00),vv5(ip00), &
         vv6(ip00),vv7(ip00),mmut(ip00)
!
    n0=npc(l)
!
    do m=1,nnc(l)
       vv1(m)=0.
       vv2(m)=0.
       vv3(m)=0.
       vv4(m)=0.
       vv5(m)=0.
       mmut(m)=0.
    enddo
    do m=1,nnc(l)
       n=m+n0
       vv1(m)=v(n,1)
       vv2(m)=v(n,2)
       vv3(m)=v(n,3)
       vv4(m)=v(n,4)
       vv5(m)=v(n,5)
    enddo
!
    if (equat(1:2).eq.'ns') then
       do m=1,nnc(l)
          n=m+n0
          mmut(m)=mut(n)
       enddo
    endif
!
    if(equat(6:7).eq.'ke') then
       do m=1,nnc(l)
          vv6(m)=0.
          vv7(m)=0.
       enddo
       do m=1,nnc(l)
          n=m+n0
          vv6(m)=v(n,6)
          vv7(m)=v(n,7)
       enddo
    endif
!
    return
  end subroutine cccca
end module mod_cccca
