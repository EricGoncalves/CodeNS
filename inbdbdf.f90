module mod_inbdbdf
  implicit none
contains
  subroutine inbdbdf(clmf,ibddim,nvbc,vbc)
!
!***********************************************************************
!
!     ACT
!_A    Adimensionnement des valeurs affectees en condition limite
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use definition
    use sortiefichier
    implicit none
    integer          :: ibddim,  nvbc
    double precision :: alf0,bet0,pdim,tdim, vbc
!
!-----------------------------------------------------------------------
!
    character(len=4 ) :: clmf
    dimension vbc(ista*lsta)
!
    if(ibddim.eq.0) then
       pdim=1.
       tdim=1.
    elseif(ibddim.eq.1) then
       pdim=pnz
       tdim=tnz
    endif
!
    if(clmf(1:3).eq.'idd') then
       vbc(1)=vbc(1)/pdim
       vbc(2)=vbc(2)/tdim
       if (nvbc.gt.2) then
          alf0=vbc(3)*atan(1.)/45.
          bet0=vbc(4)*atan(1.)/45.
       else
          alf0=0.
          bet0=0.
       endif
       vbc(3)=cos(alf0)*cos(bet0)
       vbc(4)=-cos(alf0)*sin(bet0)
       vbc(5)=sin(alf0)
       nvbc=5
!
    elseif(clmf(1:3).eq.'idi') then
       vbc(1)=vbc(1)/pdim
       vbc(2)=vbc(2)/tdim
       nvbc=2
!
    elseif(clmf.eq.'gli1') then
       nvbc=0
!
    else if (clmf.eq.'glis') then
       nvbc=0
!
    else if (clmf.eq.'pari') then
       vbc(1)=vbc(1)/tdim*ta1
       nvbc=1
!
    else if (clmf.eq.'para') then
       nvbc=0
!
    else if (clmf.eq.'prec') then
       vbc(1)=vbc(1)/pdim*pa1
       nvbc=1
!
    else if (clmf(1:4).eq.'prd ') then
       vbc(1)=vbc(1)/pdim*pa1
       nvbc=1
!
    else if (clmf(1:4).eq.'prdp') then
       vbc(1)=vbc(1)/pdim*pa1
       nvbc=1
!
    else if (clmf.eq.'extr') then
       nvbc=0
!
    else if (clmf(1:3).eq.'nrd') then
       nvbc=3
!
    else if (clmf(1:3).eq.'vrt') then
       nvbc=3
!
    else if (clmf(1:3).eq.'sym') then
       nvbc=0
!
    else if (clmf(1:3).eq.'rec') then
       nvbc=0
!
    else if (clmf(1:3).eq.'rnc') then
       nvbc=0
!
    else if (clmf(1:2).eq.'rc') then
       nvbc=0
!
    else if (clmf(1:4).eq.'rien') then
       nvbc=0
!
    else if (clmf(1:3).eq.'axe') then
       nvbc=0
!c    condition aux limites lois de paroi
!c    lp2 --> lois de paroi standard - parois adiabatiques
!c    lp3 --> lois de paroi standard - parois isothermes
!c    lois de paroi en parois adiabatiques
    elseif((clmf(1:3).eq.'lp2').or.(clmf(1:3).eq.'lp4').or. &
         (clmf(1:3).eq.'lp5')) then
       nvbc=0
!c    lois de paroi en parois isothermes
    else if (clmf(1:3).eq.'lp3') then
       vbc(1)=vbc(1)/tdim*ta1
       nvbc=1
    else if (clmf(1:4).eq.'choc') then
       nvbc=0
    else if (clmf(1:4).eq.'acou') then
       nvbc=0
    else if (clmf(1:4).eq.'debi') then
       nvbc=1
!
    else
       if (kimp.ge.1) then
          write(imp,'(a,a)') 'condition limite inconnue : ',clmf
       endif
       stop
    endif
!
    return
  end subroutine inbdbdf
end module mod_inbdbdf
