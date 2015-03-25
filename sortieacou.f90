module mod_sortieacou
  implicit none
contains
  subroutine sortieacou(l,t)
!
!***********************************************************************
!
!_DA  DATE_C :  mai 2002 -- AUTEUR : SINUMEF Eric Goncalves
!
!     ACT
!_A   Sortie pour calculs aeroacoustiques.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use chainecarac
    use maillage
    implicit none
  integer          ::    i,  i1,i1m1,  i2,i2m1
  integer          :: indc,   j,  j1,j1m1,  j2
  integer          :: j2m1,   k,  k1,k1m1,  k2
  integer          :: k2m1,   l,   m,   n, n0c
  integer          ::  nft, nid,nijd, njd
  double precision :: t(ip11,ip60)
  logical          :: ouvert
!
!-----------------------------------------------------------------------
!
    character*40 nom
    character(len=1 ) :: c
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
    nft=80
    inquire(nft,opened=ouvert,name=nom)
    if(.not. ouvert) then
!        Premier appel. Ouverture fichier et ecriture entete.
       open(nft,file='sortieacou',form='formatted')
!
       write(nft,'(''TITLE='',a1,a80,a1)')c,titrt1,c
       write(nft,'(''VARIABLES = '',a1,3(a,a1,'', '',a1),a,a1)') &
            c,'rho',c, c,'rho_u',c, c,'rho_w',c, c,'rho_E',c
       write(nft,'("ZONE F=POINT, I=",i3," J=",i3)')i2m1,j2m1
    endif
!
    k=1
    j=32
    do i=i1,i2m1
       n=indc(i,j,k)
       m=n-n0c
       write(nft,'(4(1pe15.6))') &
            t(n,1),t(n,2),t(n,4),t(n,5)
    enddo
!
    return
  end subroutine sortieacou
end module mod_sortieacou
