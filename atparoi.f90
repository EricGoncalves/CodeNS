module mod_atparoi
  implicit none
contains
  subroutine atparoi(l,nfbe)
!
!***********************************************************************
!
!     ACT
!_A    recherche du numero externe de frontiere correspondant a une paroi
!_A    dans le domaine l.
!_A    Retourne -1 s'il n'y a pas de paroi
!
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I                                        dans tab toutes front
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!
!     OUT
!_O    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_O                                        cell. interieure adjacente a la front
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use boundary
    use maillage
    implicit none
    integer          ::    l,   n,  nf,nfbe,nfbi
!
!-----------------------------------------------------------------------
!
    nf=0
    do n=1,mtbx
       nfbi=nfei(n)
       if(ndlb(nfbi).eq.l) then
!
!         la frontiere "n" est dans le domaine "l"
!
          if((cl(nfbi)(1:2).eq.'pa').or.(cl(nfbi)(1:2).eq.'lp')) then
             nfbe=n
             nf=nf+1
          end if
       end if
    enddo
    if(nf.eq.0) then
!
!       paroi non trouvee
!
       nfbe=-1
    else if(nf.ge.2) then
!
!       plus de 2 parois dans le bloc
!       A FAIRE
!
!       FAIRE message ou verifier a l'initialisation
!       --------------------------------------------
!
       write(imp,'(/,''!!!atparoi: plus de 2 paroi dans domaine '',''STOP'')')
       stop
    end if
!
    return
  end subroutine atparoi
end module mod_atparoi
