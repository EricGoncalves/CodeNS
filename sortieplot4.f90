      subroutine sortieplot(x,y,z,l,t,ps,dist,mu,mut)
!
!***********************************************************************
!
!_DA  DATE_C :  novembre 2001 -- AUTEUR : SINUMEF Eric Goncalves
!
!     ACT
!_A   Sortie tecplot pour calculs instationnaires.
!_A    Taux de vide, vitesse, pression.
!
!     VAL
!_L    titrt1     : com char             ; titre du calcul
!_L    c          :     char             ; caractere "
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use sortiefichier
      use chainecarac
      use proprieteflu
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
double precision :: x
double precision :: y
double precision :: z
integer :: l
double precision :: t
double precision :: ps
double precision :: dist
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: m
integer :: n
integer :: n0c
integer :: nid
integer :: nijd
integer :: njd
double precision :: xcc
double precision :: ycc
double precision :: zcc
!
!-----------------------------------------------------------------------
!
      logical ouvert
      character(len=1 ) :: c
      real u,v,w,e,mu,mut
!
      dimension t(ip11,ip60)
      dimension x(ip21),y(ip21),z(ip21),ps(ip11)
      dimension dist(ip12),mu(ip12),mut(ip12)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0c=npc(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!     double cote
      c=char(34)
!
        write(sec,'(''TITLE='',a1,a50,a1)')c,titrt1,c
        write(sec,'(''VARIABLES = '',a1,13(a,a1,'', '',a1),a,a1)') &
          c,'x',c, c,'y',c, c,'u',c, c,'v',c, c,'rho',c, c,'Pstat',c, &
          c,'dist',c, c,'kk',c, c,'w',c, c,'mu',c, c,'mut',c, c,'mutsmu',c
        write(sec,'("ZONE F=POINT, I=",i3," J=",i3)')j2m1,i2m1
!
        do k=k1,k2m1
!         do j=j1,j2m1
          do i=i1,i2m1
           do j=j1,j2m1
           n=indc(i,j,k)
           m=n-n0c
!
           xcc=(x(n)     +x(n     +1)+x(n     +nid)+x(n     +nid+1) &
                +x(n+nijd)+x(n+nijd+1)+x(n+nijd+nid)+x(n+nijd+nid+1))* &
                0.125
           ycc=(y(n)     +y(n     +1)+y(n     +nid)+y(n     +nid+1) &
                +y(n+nijd)+y(n+nijd+1)+y(n+nijd+nid)+y(n+nijd+nid+1))* &
                0.125
           zcc=(z(n)     +z(n     +1)+z(n     +nid)+z(n     +nid+1) &
                +z(n+nijd)+z(n+nijd+1)+z(n+nijd+nid)+z(n+nijd+nid+1))* &
                0.125
!
           u=t(n,2)/t(n,1)
           v=t(n,3)/t(n,1)
           w=t(n,4)/t(n,1)
           e=t(n,5)/t(n,1)
!
           write(sec,'(12(1pe15.6))') &
             xcc,ycc,u,v,t(n,1),ps(n),dist(n),t(n,6),t(n,7),mu(n),mut(n),mut(n)/mu(n)
          enddo
         enddo
        enddo
!
      return
      end
