module mod_defdfst
  implicit none
contains
  subroutine defdfst
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des valeurs par defaut necessaires pour l'action
!_A    dfst.
!
!***********************************************************************
!
    use para_fige
    use constantes
    use kcle
    use definition
    implicit none
    integer          :: lst,nst
!
!-----------------------------------------------------------------------
!
    do lst=1,lsta
       do nst=1,nsta
          varst(nst,lst)=reelmx
          kvarst(nst,lst)=1
       enddo
    enddo
    do nst=1,nsta
       varst(nst,1)=1.
       varst(nst,2)=1.
       varst(nst,3)=300.
       varst(nst,4)=1.
       varst(nst,5)=0.
       varst(nst,6)=0.
       varst(nst,7)=0.
    enddo
!
    return
  end subroutine defdfst
end module mod_defdfst
