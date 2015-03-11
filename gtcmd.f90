module mod_gtcmd
implicit none
contains
      subroutine gtcmd(command,lgcmd)
!
!***********************************************************************
!
!     ACT
!_A    Reconstitution d'une commande a partir
!_A    eventuellement de plusieurs lignes lues
!_A    (& etant le caractere de suite)
!_A    et determination de sa longueur.
!
!-----parameters figes--------------------------------------------------
!
      use para_fige
      use sortiefichier
      use constantes 
implicit none
integer :: lgcmd
integer :: icmd
integer :: ilin
integer :: ipos
integer :: ksuite
integer :: lglin
integer :: lgt
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: command
      character(len=132)  :: line
!
      do icmd=1,lgcmdx
       command(icmd:icmd)=' '
      enddo
      lgcmd=0
!
      lglin=0
      ksuite=0
      do while(ksuite.lt.2)
      read(lec,'(a)') line
      if(line(1:1).eq.'#') then
        lgcmd=1
        command(lgcmd:lgcmd)='#'
        return
      endif
!
      ksuite=0
      do ipos=1,linx
      if(line(ipos:ipos).ne.' ') lglin=ipos
      if(line(ipos:ipos).eq.'&') then
        ksuite=1
        exit
      endif
      enddo
!
      if(lglin.eq.0) return
      lgt=lgcmd+lglin
      if(lgt.gt.lgcmdx) stop ' commande trop longue '
      do ilin=1,lglin
       icmd=lgcmd+ilin
       command(icmd:icmd)=line(ilin:ilin)
      enddo
      lgcmd=lgcmd+lglin
      if(ksuite.eq.0) return
      if(ksuite.eq.1) command(lgcmd:lgcmd)=' '
!
      enddo
!
      return
      end subroutine
end module
