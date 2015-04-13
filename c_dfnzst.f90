module mod_c_dfnzst
  implicit none
contains
  subroutine c_dfnzst(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dfnzst.
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    gam        : com real             ; rapport des chaleurs specifiques
!
!     OUT
!_O    rgp        : com real             ; constante des gaz parfaits adim
!_O    cp         : com real             ; chal spec a pres cste adim
!_O    cv         : com real             ; chal spec a vol cst adim
!
!
!     COM
!_C   Vitesse du son et longueur de normalisation: anz,dnz
!_C   adimensionnent la vitesse de rotation.
!_C
!_C   Reynolds et temperature de normalisation: reynz,tnz
!_C   necessaires pour calcul de viscosite moleculaire, zvismo.
!_C
!_C   chaleur specifique a pression constante: cp
!_C   necessaire pour calcul des termes visqueux, zfluto.
!_C
!_C   chaleur specifique a masse volumique constante: cv
!_C   necessaire dans condition de paroi isotherme, clpari.
!_C
!_C   calcul de cp et cv dans l'adimensionnement choisi, avec
!_C   les notations:
!_C   rd pour la constante des gaz parfaits dimensionnee,
!_C   ra ou rgp dans le code, cette meme constante adimensionnee.
!_C   anz**2=gam*rd*tnz donc rgp=ra=rd/(anz**2/tnz)=1./gam
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use mod_tcmd_dfnzst
    use mod_b1_dfnzst
    use mod_dfnzst
    implicit none
    integer          :: imot(nmx),     nmot,   nonzst
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!$OMP MASTER
!
    call tcmd_dfnzst( &
         mot,imot,nmot, &
         nonzst)
!
    if (kimp.ge.1) then
       call b1_dfnzst(nonzst)
    endif
!
    call dfnzst(nonzst)
!
!$OMP END MASTER
    return
  end subroutine c_dfnzst
end module mod_c_dfnzst
