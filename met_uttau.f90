module mod_met_uttau
  implicit none
contains
  subroutine met_uttau( &
       ncbd, &
       nxn,nyn,nzn, &
       toxx,toxy,toxz,toyy,toyz,tozz, &
       s,utau)
!
!***********************************************************************
!
!     ACT
!_A    calcul du module de la vitesse de frottement sur toutes les facettes
!_A    des parois
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use sortiefichier
    implicit none
    integer          ::         i1,        i2,      i2m1,        j1,        j2
    integer          ::       j2m1,        k1,        k2,      k2m1,         l
    integer          ::         m0,       m0n,        m1,        mf,      mfac
    integer          ::      mfacn,       mfl,        mt,        n0,       n0c
    integer          :: ncbd(ip41),       nci,       ncj,       nck,     nfacf
    integer          ::        nid,      nijd,       njd
    double precision ::    nxn(ip42),   nyn(ip42),   nzn(ip42),s(ip11,ip60),     taunorm
    double precision ::        taupe,  toxx(ip12),  toxy(ip12),  toxz(ip12),  toyy(ip12)
    double precision ::   toyz(ip12),  tozz(ip12),  utau(ip42),         utx,        utxt
    double precision ::          uty,        utyt,         utz,        utzt
!
!-----------------------------------------------------------------------
!
!
    do mf=1,nbdko
!       Boucles sur les parois
!
       mfl=lbdko(mf)
       l=ndlb(mfl)
!       mfl  : numero interne frontiere (paroi)
!       l    : numero domaine
       n0=npn(l)
       n0c=npc(l)
       i1=ii1(l)
       i2=ii2(l)
       j1=jj1(l)
       j2=jj2(l)
       k1=kk1(l)
       k2=kk2(l)
       j2m1=j2-1
       k2m1=k2-1
       i2m1=i2-1
       nid = id2(l)-id1(l)+1
       njd = jd2(l)-jd1(l)+1
       nijd = nid*njd
       nci=1
       ncj=nid
       nck=nijd
!
       m0 =mpb(mfl)
       m0n=mpn(mfl)
       mt =mmb(mfl)
!       m0  : pointeur deniere facette frt. precedente
!       m0n : idem pour frontieres a normale stockee
!       mt  : nombre de facettes dans frontiere
!
!$OMP SIMD
       do m1=1,mt
!         boucle sur les facettes de la paroi
!
          mfac=mpb(mfl)+m1
          mfacn=m0n +m1
          nfacf=ncbd(mfac)
!                          mfac  : pointeur facette toutes frontieres
!                          mfacn : pointeur facette frontieres normale stockee
!                          mfacf : pointeur cellule frontiere fictive

!         calcul des composantes du vecteur contrainte T=sigma^n aux parois
!
          utx=toxx(nfacf)*nxn(mfacn)+toxy(nfacf)*nyn(mfacn)+ &
               toxz(nfacf)*nzn(mfacn)
          uty=toxy(nfacf)*nxn(mfacn)+toyy(nfacf)*nyn(mfacn)+ &
               toyz(nfacf)*nzn(mfacn)
          utz=toxz(nfacf)*nxn(mfacn)+toyz(nfacf)*nyn(mfacn)+ &
               tozz(nfacf)*nzn(mfacn)
!
!         projection du frottement sur la surface
!
          taunorm=utx*nxn(mfacn)+uty*nyn(mfacn)+utz*nzn(mfacn)
          utxt   =utx-taunorm*nxn(mfacn)
          utyt   =uty-taunorm*nyn(mfacn)
          utzt   =utz-taunorm*nzn(mfacn)
!
          taupe=sqrt(utxt**2+utyt**2+utzt**2)
          utau(mfacn)=sqrt(taupe/s(nfacf,1))
!
!         fin boucle sur les facetttes
       enddo
!       fin boucle sur les parois
    enddo
!
    return
  end subroutine met_uttau
end module mod_met_uttau
