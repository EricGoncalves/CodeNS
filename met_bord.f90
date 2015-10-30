module mod_met_bord
  implicit none
contains
  subroutine met_bord( &
       kpst, &
       bceqt, &
       mnc,ncin,mnr,xnr,ynr,znr,ncbd, &
       v,utau,mu)
!
!***********************************************************************
!
!     ACT
!      Condition aux limites pour quantites turbulentes
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use boundary
    use definition
    use sortiefichier
    use mod_met_rbve
    use mod_met_rbvc
    use mod_met_rbvr
    implicit none
    integer          ::       kpst,         m,      m0ns,        mb,       mfb
    integer          ::  mnc(ip43), mnr(ip44),        mt,        n2,        nc
    integer          :: ncbd(ip41),ncin(ip41),    nfacns,        ni,        nl
    integer          ::         no
    double precision :: bceqt(ip41,neqt),        mu(ip12),      utau(ip42),    v(ip11,ip60),       xnr(ip44)
    double precision ::        ynr(ip44),       znr(ip44)
!
!-----------------------------------------------------------------------
!
!
!
    do mfb=1,mtcx
       lbd(mfb)=nfbc(mfb)
    enddo
    nbd=mtcx
    call met_rbvc(v, ncbd,ncin,mnc) ! TODO before or after

    nbd=1
!
    do no=1,mtbx
!
       mfb=nba(no)
       lbd(1)=mfb
!
!...condition d' injection , direction de la vitesse imposee
       if (cl(mfb)(1:3).eq.'idd') then
          call met_rbve(v,ncin,ncbd)
!
!           traitement k-e impose en entree de domaine
          mt=mmb(mfb)
          if (cl(mfb)(4:4).eq.' ') then
!             valeurs d'injection egales aux valeurs a l'infini
             do m=1,mt
                mb=mpb(mfb)+m
                nc=ncbd(mb)
                v(nc,6)=rokinf
                v(nc,7)=roeinf
             enddo
          elseif (cl(mfb)(4:4).eq.'k') then
!             valeurs d'injection donnees par fichier auxiliaire

             do m=1,mt
                mb=mpb(mfb)+m
                nc=ncbd(mb)
                v(nc,6)=bceqt(mb,6)
                v(nc,7)=bceqt(mb,7)
             enddo
          else
             write(imp,'("!!!met_bord: frontiere ",a4," non prevue. Possible: idd ou iddk")')cl(mfb)
             stop
          end if
!
!...condition d' injection , (direction de la vitesse)/dt=0
       else if (cl(mfb)(1:3).eq.'idi') then
          call met_rbve(v,ncin,ncbd)
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nc=ncbd(mb)
             v(nc,6)=rokinf
             v(nc,7)=roeinf
          enddo
!
!...condition de glissement
       else if (cl(mfb).eq.'gli1') then
          call met_rbve(v,ncin,ncbd)
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...condition de glissement
       else if (cl(mfb).eq.'glis') then
          call met_rbve(v,ncin,ncbd)
!
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...condition de paroi isotherme
       elseif (cl(mfb)(1:4).eq.'pari') then
          mt=mmb(mfb)
          if (kparoi.eq.0) then
!             traitement standard
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                v(nl,6)=0.
                v(nl,7)=0.
             enddo
          elseif (kparoi.eq.1) then
!            modele k-omega de base
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                ni=ncin(mb)
                v(nl,6)=0.
                v(nl,7)=v(ni,7)
             enddo
          elseif (kparoi.eq.2) then
!            calcul de omega a la paroi
             m0ns=mpn(mfb)
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                ni=ncin(mb)
                nfacns=m0ns+m
                v(nl,6)=0.
                v(nl,7)=(v(ni,1)*utau(nfacns)*(50./rkplus))**2/mu(ni)
             enddo
          else
             write(imp,'("!!!met_bord: kparoi non prevu")')
             stop
          endif
!
!...condition de paroi adiabatique
       else if (cl(mfb)(1:4).eq.'para') then
          mt=mmb(mfb)
          if (kparoi.eq.0) then
!             traitement standard
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                v(nl,6)=0.
                v(nl,7)=0.
             enddo
          elseif (kparoi.eq.1) then
!             modele k-omega de base
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                ni=ncin(mb)
                v(nl,6)=0.
                v(nl,7)=v(ni,7)
             enddo
          elseif (kparoi.eq.2) then
!             calcul de omega a la paroi
!             modele k-omega rugueux et bas Reynolds
             rkplus=5.
             m0ns=mpn(mfb)
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                ni=ncin(mb)
                nfacns=m0ns+m
                v(nl,6)=0.
                v(nl,7)=(v(ni,1)*utau(nfacns)*(50./rkplus))**2/mu(ni)
             enddo
          else
             write(imp,'("!!!met_bord: kparoi non prevu")')
             stop
          endif
!
!...condition de pression , pression imposee et vitesse extrapolee
       else if (cl(mfb).eq.'prec') then
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...condition de pression , pression imposee(relations de compatibilite)
       else if (cl(mfb)(1:3).eq.'prd') then
          call met_rbve(v,ncin,ncbd)
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...condition de debit, debit impose (relations de compatibilite)
       else if (cl(mfb)(1:4).eq.'debi') then
          call met_rbve(v,ncin,ncbd)
!
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...condition d'extrapolation
       else if (cl(mfb).eq.'extr') then
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...condition de non reflexion ou de choc
       else if((cl(mfb)(1:3).eq.'nrd').or. &
            (cl(mfb)(1:4).eq.'choc')) then
          call met_rbve(v,ncin,ncbd)
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=rokinf
             v(nl,7)=roeinf
          enddo
!
!----condition de vorticite (2D transsonique)-
       else if (cl(mfb)(1:3).eq.'vrt') then
          call met_rbve(v,ncin,ncbd)
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=rokinf
             v(nl,7)=roeinf
          enddo
!
!...condition de symetrie par rapport aux facettes frontieres
       else if (cl(mfb)(1:3).eq.'sym') then
          mt=mmb(mfb)
          do m=1,mt
             mb=mpb(mfb)+m
             nl=ncbd(mb)
             n2=ncin(mb)
             v(nl,6)=v(n2,6)
             v(nl,7)=v(n2,7)
          enddo
!
!...raccord recouvrant
       else if (cl(mfb)(1:3).eq.'rec') then
          call met_rbvr( &
               ncbd,mnr,xnr,ynr,znr, &
               v, &
               ncin)
!
!...raccord non coincident
       else if (cl(mfb)(1:3).eq.'rnc') then
          call met_rbvr( &
               ncbd,mnr,xnr,ynr,znr, &
               v, &
               ncin)
!
!...raccord coincident
       else if (cl(mfb)(1:2).eq.'rc') then
!          call met_rbvc(v, ncbd,ncin,mnc)
!
!...affectation de valeurs extrapolees en vue de post-traitement
       else if (((cl(mfb).eq.'rien').or.(cl(mfb).eq.'axe ')).and. &
            (kpst.eq.1)) then
          call met_rbve(v,ncin,ncbd)
!
!...lois de paroi -
       elseif (cl(mfb)(1:2).eq.'lp') then
          mt=mmb(mfb)
          if(kparoi.eq.0) then
!            traitement standard
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                v(nl,6)=0.
                v(nl,7)=0.
             enddo
          elseif ((kparoi.eq.1).or.(kparoi.eq.2)) then
!           modeles k-omega et k-omega bas Reynolds
             do m=1,mt
                mb=mpb(mfb)+m
                nl=ncbd(mb)
                ni=ncin(mb)
                v(nl,6)=0.
                v(nl,7)=v(ni,7)
             enddo
          else
             write(imp,'("!!!met_bord: kparoi non prevu")')
             stop
          endif
!
!...aucun traitement n'est effectue
       else if (((cl(mfb).eq.'rien').or.(cl(mfb).eq.'axe ')).and. &
            (kpst.ne.1)) then
!
       else
          if (kimp.ge.1) then
             write(imp,'(a,a)') 'condition limite inconnue : ',cl(mfb)
          endif
          stop
       endif
!
    enddo

    do mfb=1,mtcx
       lbd(mfb)=nfbc(mfb)
    enddo
    nbd=mtcx
    call met_rbvc(v, ncbd,ncin,mnc) ! TODO before or after
!
    return
  end subroutine met_bord
end module mod_met_bord
