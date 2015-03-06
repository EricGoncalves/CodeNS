      subroutine b1_dfnm
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfnm'.
!
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    nfi        : com int              ; nbr de rangees de pts fictifs
!
!-----------------------------------------------------------------------
!
      use sortiefichier
      use maillage
      use schemanum
      use kcle
      use proprieteflu
!
!-----------------------------------------------------------------------
!
      character *1316 form
      character *24 cnfi, &
                    ckfmg,clgx,ckcg, &
                    ckdualns,ctol,ctolke,cniter,cnitur, &
                    cischema,cmuscl,cilim,cxk,ckprec,ccte,ckvisq
!
      call convich(knfi,cnfi)
      call convich(kkfmg,ckfmg)
      call convich(klgx,clgx)
      call convich(kkcg,ckcg)
      call convich(kischema,cischema)
      call convich(kmuscl,cmuscl)
      call convich(kilim,cilim)
      call convich(kxk,cxk)
      call convich(kkprec,ckprec)
      call convich(kcte,ccte)
      call convich(kkvisq,ckvisq)
      call convich(kkdualns,ckdualns)
      call convich(ktol,ctol)
      call convich(kniter,cniter)
      call convich(ktolke,ctolke)
      call convich(knitur,cnitur)
!
           form='(/,2x,''definition du schema numerique'',/' &
                 //'2x,''------------------------------'',/' &
                 //'2x,''nb de pts fictifs [nfi]  : '',11x,i5,2x,a)'
      write(imp,form) nfi,cnfi
!
            form='(2x,''full multigrid [kfmg]    : '',11x,i5,2x,a/' &
                //'2x,''               [lgx ]    : '',11x,i5,2x,a/' &
                //'2x,''               [kcg ]    : '',11x,i5,2x,a)'
      write(imp,form) kfmg,ckfmg,lgx,clgx,kcg,ckcg
!
            form='(2x,''type de schema [ischema] : '',11x,i5,2x,a)'
      write(imp,form) ischema,cischema
!
         form='(2x,''extrapolation MUSCL [muscl]  : '',11x,i5,2x,a/' &
             //'2x,''                    [ilim ]  : '',11x,i5,2x,a/' &
             //'2x,''                    [xk   ]  : '',4x,e12.6,2x,a)'
      write(imp,form) muscl,cmuscl,ilim,cilim,xk,cxk
!
         form='(2x,''preconditionnement [kprec]  : '',11x,i5,2x,a/' &
             //'2x,''                   [cte  ]  : '',4x,e12.6,2x,a/' &
             //'2x,''                   [kvisq]  : '',11x,i5,2x,a)'
      write(imp,form) kprec,ckprec,cte,ccte,kvisq,ckvisq
!
      if(kfmg.eq.3) then
            form='(2x,''temps dual   [kdualns ]  : '',11x,i5,2x,a/' &
                //'2x,''             [tol     ]  : '',4x,e12.6,2x,a/' &
                //'2x,''             [niter   ]  : '',4x,i5,2x,a/' &
                //'2x,''             [tolke   ]  : '',4x,e12.6,2x,a/' &
                //'2x,''             [nitur   ]  : '',4x,i5,2x,a)'
       write(imp,form) kdualns,ckdualns,tol,ctol,niter,cniter, &
            tolke,ctolke,nitur,cnitur
      endif
!
      return
      end
