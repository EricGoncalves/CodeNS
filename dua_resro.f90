module mod_dua_resro
  implicit none
contains
  subroutine dua_resro( &
       icyc,ncyc,img, &
       u0,v,dt)
!
!--------------------------------------------------------------------
!
!_D   DATE: novembre 2001 - Eric Goncalves / SINUMEF
!
!_A   ACT: Calcul du residu (norme L2) pour la densite
!
!_O   OUT: durmy2   ; residus quadratiques moyens de rho sur le domaine
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use schemanum
    use sortiefichier
    use mod_mpi,only:rank,sum_mpi
    implicit none
    integer          ::    i,  i1,  i2,i2m1,icyc
    integer          ::  img,   j,  j1,  j2,j2m1
    integer          ::    k,  k1,  k2,k2m1,   l
    integer          ::   lm,   n,  n0,ncyc, nid
    integer          :: nijd, njd
    double precision ::      dt(ip11),       durmy2,         resr,u0(ip11,ip60), v(ip11,ip60)
!
!--------------------------------------------------------------------
!
!
    durmy2 = 0.
    do l=1,lzx
       lm=l+(img-1)*lz
       n0=npc(lm)
       i1=ii1(lm)
       i2=ii2(lm)
       j1=jj1(lm)
       j2=jj2(lm)
       k1=kk1(lm)
       k2=kk2(lm)
       i2m1=i2-1
       j2m1=j2-1
       k2m1=k2-1
!
       nid=id2(lm)-id1(lm)+1
       njd=jd2(lm)-jd1(lm)+1
       nijd=nid*njd
!
       do k=k1,k2m1
          do j=j1,j2m1
             do i=i1,i2m1
                n=ind(i,j,k)
                resr=(v(n,1)-u0(n,1))/dt(n)
                durmy2=durmy2+resr*resr
             enddo
          enddo
       enddo
    enddo        ! Fin boucle domaines - Grille fine
    call sum_mpi(durmy2)
!
      if(ncyc.eq.1) then
       resno1=sqrt(durmy2)
       if(abs(resno1)<tiny(1.)) resno1=1.
      endif
      resite=sqrt(durmy2)/resno1
!     resite=sqrt(durmy2)
    if (rank==0) then
      open(sor3 ,file='resro',position="append")
      write(sor3,'(1x,i6,1x,i6,1x,e13.6)') ncyc,icyc,resite
      close(sor3)
    endif
!
    return
  contains
    function    ind(i,j,k)
      implicit none
      integer          ::   i,ind,  j,  k
      ind=n0+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function ind
  end subroutine dua_resro
end module mod_dua_resro
