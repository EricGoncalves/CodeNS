module mod_idirch
implicit none
contains
      function idirch(ch)
implicit none
integer :: idir
!
!***********************************************************************
!
!     ACT
!_A    Transformation d'une chaine de caractere i1..k2 en entier.
!_A      i1 --> +1     i2--> -1
!_A      j1 --> +2     j2--> -2
!_A      k1 --> +3     k2--> -3
!
!     INP
!_I    ch         : arg char             ; type de plan
!
!     OUT
!_O    valeur de la fonction
!
!-----------------------------------------------------------------------
!
      character(len=2 ) :: ch
      integer idirch
!
      idirch = 0
!
      if(ch(2:2).eq.'2') then
       idir = -1
      elseif(ch(2:2).eq.'1') then
       idir = +1
      else
       write(*,1000)
       stop
      end if
!
      if(ch(1:1).eq.'i') then
       idirch = 1*idir
      elseif(ch(1:1).eq.'j') then
       idirch = 2*idir
      elseif (ch(1:1).eq.'k') then
       idirch = 3*idir
      else
       write(*,1000)
       stop
      end if
!
      if (idirch.eq.0) then
       write(*,1000)
       stop
      end if
!
 1000 format(10x,'caracter - type of subdomain face incomprehensible')
!
      return
      end function
end module
