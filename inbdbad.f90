module mod_inbdbad
  implicit none
contains
  subroutine inbdbad( &
       mfl,nv,adim, &
       chvar)
!
!***********************************************************************
!
!     ACT
!_A     calcul des correspondances conditions limites <-> numero
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use definition
    use boundary
    implicit none
    integer          :: mfl, nv
    double precision :: adim
!
!-----------------------------------------------------------------------
!
    character(len=20) ::  chvar
!
    nv=1
    adim=1.
!
    if (cl(mfl)(1:3).eq.'idd' .or. &
         cl(mfl)(1:3).eq.'idi' .or. &
         cl(mfl)(1:4).eq.'dist') then
!
       if (chvar(1:2).eq.'PI' .or.chvar(1:2).eq.'pi' ) then
          nv=1
       endif
       if (chvar(1:2).eq.'TI' .or.chvar(1:2).eq.'ti' ) then
          nv=2
       endif
       if (chvar(1:3).eq.'D0X'.or.chvar(1:3).eq.'d0x') then
          nv=3
       endif
       if (chvar(1:3).eq.'D0Y'.or.chvar(1:3).eq.'d0y') then
          nv=4
       endif
       if (chvar(1:3).eq.'D0Z'.or.chvar(1:3).eq.'d0z') then
          nv=5
       endif
       if (chvar(1:3).eq.'ROK'.or.chvar(1:3).eq.'rok') then
          nv=6
       endif
       if (chvar(1:4).eq.'ROEP'.or.chvar(1:4).eq.'roep') then
          nv=7
       endif
!
!
    else if (cl(mfl)(1:4).eq.'gli1'.or. &
         cl(mfl)(1:4).eq.'glis') then
!
    else if (cl(mfl)(1:4).eq.'prec'.or. &
         cl(mfl)(1:4).eq.'prd '.or. &
         cl(mfl)(1:4).eq.'prdp') then
!
       if (chvar(1:2).eq.'PS' .or.chvar(1:2).eq.'ps' ) then
          nv=1
          adim=pnz/pa1
       endif
!
    else if (cl(mfl)(1:4).eq.'pari') then
!
       if (chvar(1:2).eq.'TP' .or.chvar(1:2).eq.'tp' ) then
          nv=1
          adim=tnz/ta1
       endif
!
    else if (cl(mfl)(1:4).eq.'para') then
!
    else if (cl(mfl)(1:4).eq.'extr') then
!
    else if (cl(mfl)(1:3).eq.'nrd' ) then
!
    else if (cl(mfl)(1:3).eq.'vrt' ) then
!
    else if (cl(mfl)(1:3).eq.'sym' ) then
!
    else if (cl(mfl)(1:3).eq.'rec' ) then
!
    else if (cl(mfl)(1:3).eq.'rnc' ) then
!
    else if (cl(mfl)(1:2).eq.'rc'  ) then
!
    else if (cl(mfl)(1:4).eq.'rien') then
!
    else if (cl(mfl)(1:3).eq.'axe' ) then
!     condition aux limites lois de paroi
!     lp2 --> lois de paroi standard - parois adiabatiques
!     lp3 --> lois de paroi standard - parois isothermes
!     lois de paroi en parois adiabatiques
    elseif((cl(mfl)(1:3).eq.'lp2').or.(cl(mfl)(1:3).eq.'lp4').or. &
         (cl(mfl)(1:3).eq.'lp5')) then
    else if(cl(mfl)(1:3).eq.'lp3') then
       if((chvar(1:2).eq.'TP').or.(chvar(1:2).eq.'tp')) then
          nv=1
          adim=tnz/ta1
       endif
    else if(cl(mfl)(1:4).eq.'choc') then
    else if(cl(mfl)(1:4).eq.'acou') then
    else if(cl(mfl)(1:4).eq.'debi') then
       if (chvar(1:2).eq.'QM' .or.chvar(1:2).eq.'qm' ) then
          nv=1
       endif
!
    endif
!
    return
  end subroutine inbdbad
end module mod_inbdbad
