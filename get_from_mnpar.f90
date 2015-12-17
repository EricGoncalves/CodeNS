module mod_get_from_mnpar
  implicit none
contains
  subroutine get_from_mnpar(var_out,var_in,mnpar,mnpar2)
    use mod_mpi
    use para_var
    use boundary
    use tools
    implicit none
    double precision,dimension(ip42),intent(in)  :: var_in
    double precision,dimension(ip12),intent(out) :: var_out
    integer,dimension(ip12)         ,intent(in)  :: mnpar,mnpar2
    
    integer         ,allocatable :: old_ranks(:) 
    double precision,allocatable :: buff(:)
    integer :: bcg,proc,new_proc,new_comm,new_rank,new_NPROCS,mbmx,m0b,mb,mbb,nc,IERR
    logical :: test,test2


    do bcg=1,num_bcg
      test=sum(mnpar2,mask=(mnpar2==bcg))/=0 ! at least one of my points is concerned
      call LOR_MPI(test,test2)
      if (test2) then
        proc=bcg_to_proc(bcg) ! I'm the owner

        call barrier ! TODO, usefull ?

        !prepare new communicator with only concerned procs
        
        if (test.or.rank==proc) then
          call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1, rank, NEW_COMM, IERR)
          CALL MPI_COMM_RANK(NEW_COMM, new_RANK, IERR)
          CALL MPI_COMM_SIZE(NEW_COMM, new_NPROCS, IERR)
          call reallocate(old_ranks,new_NPROCS)
          call gather(rank,old_ranks,comm=new_comm)
          
          new_proc=minloc(old_ranks, 1, mask=(old_ranks==proc))-1

          !prepare buffer
          mbmx=0
          if(rank==proc) mbmx=mmb(bcg_to_bcl(bcg))
          call bcast(mbmx,new_proc,comm=new_comm)
          call reallocate(buff,mbmx)
          
          !fill buffer
          if(rank==proc) then
             m0b=mpn(bcg_to_bcl(bcg))
             do mb=1,mbmx
                mbb=m0b+mb
                buff(mb)=var_in(mbb)
             enddo
          endif
          
          !send buffer
          call bcast(buff,new_proc,comm=new_comm)

          !fill var_out
          
          do nc = 1,ip12
             if(mnpar2(nc)==bcg) var_out(nc)= buff(mnpar(nc))
          enddo
          call barrier ! TODO, usefull ?
          call MPI_COMM_FREE(NEW_COMM, IERR)
        else
          call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, rank, NEW_COMM, IERR)
          call barrier ! TODO, usefull ?
          call MPI_COMM_FREE(NEW_COMM, IERR)
        endif
      endif
   enddo
  end subroutine get_from_mnpar
end module mod_get_from_mnpar
