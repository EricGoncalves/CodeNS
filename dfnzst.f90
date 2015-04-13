module mod_dfnzst
  implicit none
contains
  subroutine dfnzst(nonzst)
!
!***********************************************************************
!
!_DA  DATE_C : 16 mars 2006 - Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul des valeurs adimensionnees de caracteristiques du gaz pour
!_A    completer l'adimensionnement.
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
!_C   rgp dans le code, cette meme constante adimensionnee.
!_C   anz**2=gam*rd*tnz donc rgp=rd/(anz**2/tnz)=1./gam
!
!_C   pour Stiffened gas, rd=(gam-1)*cp/gam
!
!***********************************************************************
!
    use para_fige
    use definition
    use kcle
    use proprieteflu
    implicit none
    integer          :: nonzst
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
    ronz=varst(nonzst,1)
    tnz =varst(nonzst,3)
    dnz =varst(nonzst,4)
    rnz=rd
    anz=sqrt(gam*rnz*tnz)
    pnz=ronz*anz**2/gam
!       anz =varst(nonzst,2)
!       pnz =ronz/gam*anz**2
!       rnz=anz**2/(gam*tnz)
!
    rgp=1./gam
    cp=1./(gam-1.)
    cv=cp/gam
!
!$OMP END MASTER
    return
  end subroutine dfnzst
end module mod_dfnzst
