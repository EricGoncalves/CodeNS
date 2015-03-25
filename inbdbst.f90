module mod_inbdbst
  implicit none
contains
  subroutine inbdbst(clmf,ibdcst,nvbc,vbc)
!
!***********************************************************************
!
!     ACT
!_A    Utilisation d'un etat thermodynamique pour appliquer les
!_A    conditions aux limites
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use proprieteflu
    use definition
    implicit none
  integer          :: ibdcst,  nvbc
  double precision ::           aist,          alst,          btst,          mast,          pist
  double precision ::           pspi,         roist,          tist,          tsti,vbc(ista*lsta)
!
!-----------------------------------------------------------------------
!
    character(len=4 ) :: clmf
!
    roist =varst(ibdcst,1)
    aist  =varst(ibdcst,2)
    tist  =varst(ibdcst,3)
    mast  =varst(ibdcst,5)
    alst  =varst(ibdcst,6)
    btst  =varst(ibdcst,7)
!      rgp   =aist**2/(gam*tist)
    rgp   =1./gam
    pist  =roist*rgp*tist
!
    if(clmf(1:3).eq.'idd') then
       vbc(1)=pist
       vbc(2)=tist
       vbc(3)=alst
       vbc(4)=btst
       nvbc=4
!
    else if (clmf(1:3).eq.'idi') then
       vbc(1)=pist
       vbc(2)=tist
       nvbc=2
!
    else if (clmf.eq.'gli1') then
       nvbc=0
!
    else if (clmf.eq.'glis') then
       nvbc=0
!
    else if (clmf.eq.'pari') then
       tsti=(1.+.5*(gam-1.)*mast**2)**(-1.)
       vbc(1)=tsti*tist
       nvbc=1
!
    else if (clmf.eq.'para') then
       nvbc=0
!
    else if (clmf.eq.'prec') then
       pspi=(1.+.5*(gam-1.)*mast**2)**(gam/(1.-gam))
       vbc(1)=pspi*pist
       nvbc=1
!
    else if (clmf(1:4).eq.'prd ') then
       pspi=(1.+.5*(gam-1.)*mast**2)**(gam/(1.-gam))
       vbc(1)=pspi*pist
       nvbc=1
!
    else if (clmf(1:4).eq.'prdp') then
       pspi=(1.+.5*(gam-1.)*mast**2)**(gam/(1.-gam))
       vbc(1)=pspi*pist
       nvbc=1
!
    else if (clmf.eq.'extr') then
       nvbc=0
!
    else if (clmf(1:3).eq.'nrd') then
       vbc(1)=mast
       vbc(2)=alst
       vbc(3)=btst
       nvbc=3
!
    else if (clmf(1:3).eq.'vrt') then
       vbc(1)=mast
       vbc(2)=alst
       vbc(3)=btst
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
!     lois de paroi en parois adiabatiques
    elseif((clmf(1:3).eq.'lp2').or.(clmf(1:3).eq.'lp4').or. &
         (clmf(1:3).eq.'lp5')) then
       nvbc=0
!     lois de paroi en parois isothermes
    elseif (clmf(1:3).eq.'lp3') then
       tsti=(1.+.5*(gam-1.)*mast**2)**(-1.)
       vbc(1)=tsti*tist
       nvbc=1
    else if (clmf(1:4).eq.'choc') then
       nvbc=0
    else if (clmf(1:4).eq.'acou') then
       nvbc=0
    else if (clmf(1:4).eq.'debi') then
       vbc(1)=1
       nvbc=1
!
    else
       if (kimp.ge.1) then
          write(imp,'(a,a)') 'condition limite inconnue : ',clmf
       endif
       stop
    end if
!
    return
  end subroutine inbdbst
end module mod_inbdbst
