module mod_atindnor
  implicit none
contains
  subroutine atindnor( &
       nfbe, &
       m10,m20,m30, &
       m1min,m1max,m2min,m2max,m3min,m3max, &
       dm1,dm2,dm3,isens3)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des donnees necessaires a la determination de la distance
!_A    a la paroi pour utilisation de la longueur de melange de Michel
!
!      PRINCIPE : Les directions "1" et "2" sont dans le plan de la frontiere,
!      --------   la direction "3" est suivant la normale. Le plan frontiere
!                 est decrit par une double boucle en "m1" et "m2" croissants.
!                 La normale est decrite depuis la frontiere jusqu'a l'autre
!                 extremite du domaine par une boucle allant de "i3min" a
!                 "i3max" par accroissement de "isens3" egal a +/-1.
!
!     INP
!_I    nfbi       : arg int              ; numero interne de frontiere
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    kd2        : com int (lt        ) ; indice max en k fictif
!
!     OUT
!_O    m10        : arg int              ; ind orig dir 1 de la frontiere
!_O    m20        : arg int              ; ind orig dir 2 de la frontiere
!_O    m30        : arg int              ; ind orig dir normale a la frontiere
!_O    m1min      : arg int              ; ind min dir 1 de la frontiere
!_O    m1max      : arg int              ; ind max dir 1 de la frontiere
!_O    m2min      : arg int              ; ind min dir 2 de la frontiere
!_O    m2max      : arg int              ; ind max dir 2 de la frontiere
!_O    m3min      : arg int              ; ind min dir 3 de la frontiere
!_O    m3max      : arg int              ; ind max dir 3 de la frontiere
!_O    dm1        : arg int              ; saut entre 2 points dir 1
!_O    dm2        : arg int              ; saut entre 2 points dir 2
!_O    dm3        : arg int              ; saut entre 2 points dir 3
!_O    isens3     : arg int              ; sens var indice sur normale
!
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use maillage
    use boundary
    implicit none
    integer          ::    dm1,   dm2,   dm3,  imax,  imin
    integer          :: isens3,  jmax,  jmin,  kmax,  kmin
    integer          ::      l,   m10, m1max, m1min,   m20
    integer          ::  m2max, m2min,   m30, m3max, m3min
    integer          ::   mfbi,  nfbe,   nid,   njd
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: indmf
!
!     ..................................................................
!     l           : numero de domaine
!     nfbi        : numero interne de frontiere
!     indmf       : type indice plan frontiere
!     imin        : indice min en i de la frontiere
!     imax        : indice max en i de la frontiere
!     ..................................................................
!
    mfbi=nfei(nfbe)
    l   =ndlb(mfbi)
    indmf=indfl(mfbi)
    imin=iminb(mfbi)
    imax=imaxb(mfbi)
    jmin=jminb(mfbi)
    jmax=jmaxb(mfbi)
    kmin=kminb(mfbi)
    kmax=kmaxb(mfbi)
    nid=id2(l)-id1(l)+1
    njd=jd2(l)-jd1(l)+1
!
    if(indmf(1:1).eq.'i') then
!
!       normale perpendiculaire au plan "i=imin" ou "i=imax"
!       On a la correspondance :
!       ind=n0+(i -id1(l))*1  +(j -jd1(l))*nid+(k -kd1(l))*nij
!       ind=n0+(m3-m30   )*dm3+(m1-m10   )*dm1+(m2-m20   )*dm2
!
       m1min=jmin
       m1max=jmax
       m10  =jd1(l)
       dm1  =nid
       m2min=kmin
       m2max=kmax
       m20  =kd1(l)
       dm2  =nid*njd
!
       if(indmf.eq.'i1') then
!
!         normale dirigee de "imin" vers "imax"
!
          isens3=1
          m3min =imin
          m3max =ii2(l)
          m30   =id1(l)
          dm3   =isens3*1
       else
!
!         normale dirigee de "imax" vers "imin"
!
          isens3=-1
          m3min =imax
          m3max =ii1(l)
!         m30   =id2(l)
          m30   =id1(l)
!         dm3   =isens3*1
          dm3   =1
       end if
    else if(indmf(1:1).eq.'j') then
!
!       normale perpendiculaire au plan "j=jmin" ou "j=jmax"
!       On a la correspondance :
!       ind=n0+(i -id1(l))*1  +(j -jd1(l))*nid+(k -kd1(l))*nij
!       ind=n0+(m1-m10   )*dm1+(m2-m20   )*dm2+(m3-m30   )*dm3
!
       m1min=imin
       m1max=imax
       m10  =id1(l)
       dm1  =1
       m2min=kmin
       m2max=kmax
       m20  =kd1(l)
       dm2  =nid*njd
!
       if(indmf.eq.'j1') then
!
!         normale dirigee de "jmin" vers "jmax"
!
          isens3=1
          m3min =jmin
          m3max =jj2(l)
          m30   =jd1(l)
          dm3   =isens3*nid
       else
!
!         normale dirigee de "jmax" vers "jmin"
!
          isens3=-1
          m3min =jmax
          m3max =jj1(l)
!         m30   =jd2(l)
          m30   =jd1(l)
!         dm3   =isens3*nid
          dm3   =nid
       end if
    else if(indmf(1:1).eq.'k') then
!
!       normale perpendiculaire au plan "k=kmin" ou "k=kmax"
!       On a la correspondance :
!       ind=n0+(i -id1(l))*1  +(j -jd1(l))*nid+(k -kd1(l))*nij
!       ind=n0+(m1-m10   )*dm1+(m2-m20   )*dm2+(m3-m30   )*dm3
!
       m1min=imin
       m1max=imax
       m10  =id1(l)
       dm1  =1
       m2min=jmin
       m2max=jmax
       m20  =jd1(l)
       dm2  =nid
!
       if(indmf.eq.'k1') then
!
!         normale dirigee de "kmin" vers "kmax"
!
          isens3=1
          m3min =kmin
          m3max =kk2(l)
          m30   =kd1(l)
          dm3   =isens3*nid*njd
       else
!
!         normale dirigee de "kmax" vers "kmin"
!
          isens3=-1
          m3min =kmax
          m3max =kk1(l)
!         m30   =kd2(l)
          m30   =kd1(l)
!         dm3   =isens3*nid*njd
          dm3   =nid*njd
       end if
    end if
!
    return
  end subroutine atindnor
end module mod_atindnor
