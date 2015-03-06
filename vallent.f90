      subroutine vallent(mot,imot,nm,lent,lentd,lx,klx)
!
!***********************************************************************
!
!     ACT
!_A    Affectation des valeurs a une liste d'entiers.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use chainecarac
!
!-----------------------------------------------------------------------
!
      character *32 comment
      character *32 mot(nmx)
!
      dimension imot(nmx)
      dimension lent(nobj)
!
      do icmt=1,32
      comment(icmt:icmt)=' '
      enddo
!
        call lentier(mot(nm),imot(nm),lent,lentd,lx,klx,kerr)
        if(kerr.eq.-1) then
          comment=cd
          call synterr(mot,imot,nm,comment)
        else if(kerr.eq.0)then
          comment=cm
          call synterr(mot,imot,nm,comment)
        endif
!
      return
      end
