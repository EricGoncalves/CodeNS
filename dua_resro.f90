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
implicit none
integer :: ind
double precision :: u0
double precision :: v
double precision :: dt
integer :: i
integer :: j
integer :: k
integer :: i1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: l
integer :: lm
integer :: n
integer :: n0
integer :: nid
integer :: nijd
integer :: njd
double precision :: resr
!
!--------------------------------------------------------------------
!
      integer icyc,ncyc,img
      real durmy2
      dimension u0(ip11,ip60),v(ip11,ip60),dt(ip11)
!
      ind(i,j,k)=n0+1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
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
!
!      if(ncyc.eq.1) then
!       resno1=sqrt(durmy2)
!       if(resno1.eq.0.) resno1=1.
!      endif
!      resite=sqrt(durmy2)/resno1
      resite=sqrt(durmy2)
      write(sor3,'(1x,i6,1x,i6,1x,e13.6)') ncyc,icyc,resite
!
      return
      end subroutine
end module
