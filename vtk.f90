module mod_vtk
  implicit none
  character(len=80),allocatable :: collection(:)
  character(len=80) :: collect_file,collect_dir=""
  integer :: ncollect
  logical :: collect
  INTERFACE vtk_writer
    MODULE PROCEDURE vtk_writer_I0,vtk_writer_I1,vtk_writer_R0,vtk_writer_R1
  END INTERFACE vtk_writer
contains


  subroutine vtk_start_collection(fich,dir)
      use tools
      use mod_mpi
      implicit none
      character(*),intent(in) :: fich,dir
      
      call reallocate(collection,1)
      
      collect_file=dir//"/"//fich
      collect_dir=dir
      
      if(rank==0) then  ! write header of pvd file
        open(42,file=collect_file,status="replace")
        write(42,'(A)') '<?xml version="1.0"?>'
        write(42,'(A)') ' <VTKFile type="Collection" version="0.1">'
        write(42,'(A)') '   <Collection>'
        close(42)
      endif
      
      collect=.true.
      ncollect=0
  end subroutine vtk_start_collection

  subroutine vtk_end_collection
      use tools
      use mod_mpi
      implicit none
      integer :: i,ncollect_tot
      character(len=80),allocatable :: collection_full(:)
      
      call sum_mpi(ncollect,ncollect_tot)
      allocate(collection_full(ncollect_tot))
      if(.not.allocated(collection)) allocate(collection(1))
      call gather(collection,collection_full,ncollect)

      if(rank==0) then  ! write the rest of the pvd file
        open(42,file=collect_file,position="append")
        do i=1,ncollect_tot
          write(42,'(A,I0.3,A,I0.2,A)') '     <DataSet part="',i-1,'" file="'//trim(collection_full(i))//'"/>'
        enddo
        write(42,'(A)') '   </Collection>'
        write(42,'(A)') ' </VTKFile>'
        close(42)
      endif
      
      collect_dir=""
      collect=.false.
  end subroutine vtk_end_collection

  subroutine vtk_open(vtk_file,x,y,z,l,oi1,oi2,oj1,oj2,ok1,ok2)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      double precision,intent(in) :: x(:),y(:),z(:)
      integer,intent(in) ::l
      integer,intent(in),optional :: oi1,oi2,oj1,oj2,ok1,ok2
      character(*),intent(in) :: vtk_file

      integer :: i1,i2,j1,j2,k1,k2
      integer :: nid,njd,nijd,xyz,i,j,k
      character(len=80) :: file

      if (collect) then
        ncollect=ncollect+1
        call reallocate_s(collection,ncollect)
        collection(ncollect)=vtk_file
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif
      
      if (present(oi1)) then
        i1=oi1 ; i2=oi2        ! if given, write only a portion of the domain
        j1=oj1 ; j2=oj2
        k1=ok1 ; k2=ok2
      else
        i1=ii1(l) ; i2=ii2(l)  ! if not, write all domain
        j1=jj1(l) ; j2=jj2(l)
        k1=kk1(l) ; k2=kk2(l)
      endif
      
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd = nid*njd

      open(42,file=file,status="replace")

      ! write header

      write(42,'(A)') '<?xml version="1.0"?>'
      write(42,'(A)') '<VTKFile type="StructuredGrid"  version="0.1" byte_order="LittleEndian">'
      write(42,'(A,6I4,A)') '<StructuredGrid WholeExtent="',i1,i2,j1,j2,k1,k2,'">'
      write(42,'(A,6I4,A)') '<Piece Extent="',i1,i2,j1,j2,k1,k2,'">'

      ! write mesh

      write(42,'(A)') "<Points>"
      write(42,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      do k=k1,k2
        do j=j1,j2
         do i=i1,i2

            xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

            write(42,'(3E16.8)') x(xyz),y(xyz),z(xyz)

          enddo
        enddo
      enddo

      write(42,'(A)') '</DataArray>'
      write(42,'(A)') '</Points>'
      write(42,'(A)') '<CellData>'
      close(42)

  end subroutine vtk_open


  subroutine vtk_close(vtk_file)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      character(*),intent(in) :: vtk_file
      character(len=80) :: file

      if (collect) then
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif

      open(42,file=file,position="append")

      ! write footer

      write(42,'(A)') "</CellData>"
      write(42,'(A)') '</Piece>'
      write(42,'(A)') '</StructuredGrid>'
      write(42,'(A)') '</VTKFile>'

      close(42)

  end subroutine vtk_close




  subroutine vtk_writer_I1(vtk_file,field,name_field,l,oi1,oi2,oj1,oj2,ok1,ok2)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      integer,intent(in) :: field(:)
      integer,intent(in) ::l
      integer,intent(in),optional :: oi1,oi2,oj1,oj2,ok1,ok2
      character(*),intent(in) :: vtk_file,name_field

      integer :: i1,i2,j1,j2,k1,k2
      integer :: nid,njd,nijd,xyz,i,j,k
      character(len=80) :: vartype,file

      if (collect) then
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif

      if (present(oi1)) then
        i1=oi1 ; i2=oi2        ! if given, write only a portion of the domain
        j1=oj1 ; j2=oj2
        k1=ok1 ; k2=ok2
      else
        i1=ii1(l) ; i2=ii2(l)  ! if not, write all domain
        j1=jj1(l) ; j2=jj2(l)
        k1=kk1(l) ; k2=kk2(l)
      endif

      vartype="Int32"

      open(42,file=file,position="append")

      ! write var

      write(42,'(5A)') '<DataArray type="',trim(vartype),'" Name="',name_field,'"  NumberOfComponents="1" format="ascii">'

      do k=k1,k2-1
        do j=j1,j2-1
          do i=i1,i2-1
            nid = id2(l)-id1(l)+1
            njd = jd2(l)-jd1(l)+1
            nijd = nid*njd

            xyz =npc(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

            write(42,'(I5)') field(xyz)
          enddo
        enddo
      enddo

      ! write footer

      write(42,'(A)') '</DataArray>'

      close(42)

  end subroutine vtk_writer_I1

  subroutine vtk_writer_I0(vtk_file,field,name_field,l,oi1,oi2,oj1,oj2,ok1,ok2)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      integer,intent(in) :: field
      integer,intent(in) ::l
      integer,intent(in),optional :: oi1,oi2,oj1,oj2,ok1,ok2
      character(*),intent(in) :: vtk_file,name_field

      integer :: i1,i2,j1,j2,k1,k2
      integer :: nid,njd,nijd,xyz,i,j,k
      character(len=80) :: vartype,file

      if (collect) then
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif
      
      if (present(oi1)) then
        i1=oi1 ; i2=oi2        ! if given, write only a portion of the domain
        j1=oj1 ; j2=oj2
        k1=ok1 ; k2=ok2
      else
        i1=ii1(l) ; i2=ii2(l)  ! if not, write all domain
        j1=jj1(l) ; j2=jj2(l)
        k1=kk1(l) ; k2=kk2(l)
      endif

      vartype="Int32"

      open(42,file=file,position="append")

      ! write var

      write(42,'(5A)') '<DataArray type="',trim(vartype),'" Name="',name_field,'"  NumberOfComponents="1" format="ascii">'

      do k=k1,k2-1
        do j=j1,j2-1
          do i=i1,i2-1
            nid = id2(l)-id1(l)+1
            njd = jd2(l)-jd1(l)+1
            nijd = nid*njd

            xyz =npc(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

            write(42,'(I5)') field
          enddo
        enddo
      enddo

      ! write footer

      write(42,'(A)') '</DataArray>'

      close(42)

  end subroutine vtk_writer_I0

  subroutine vtk_writer_R1(vtk_file,field,name_field,l,oi1,oi2,oj1,oj2,ok1,ok2)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      double precision,intent(in) :: field(:)
      integer,intent(in) ::l
      integer,intent(in),optional :: oi1,oi2,oj1,oj2,ok1,ok2
      character(*),intent(in) :: vtk_file,name_field

      integer :: i1,i2,j1,j2,k1,k2
      integer :: nid,njd,nijd,xyz,i,j,k
      character(len=80) :: vartype,file

      if (collect) then
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif
      
      if (present(oi1)) then
        i1=oi1 ; i2=oi2        ! if given, write only a portion of the domain
        j1=oj1 ; j2=oj2
        k1=ok1 ; k2=ok2
      else
        i1=ii1(l) ; i2=ii2(l)  ! if not, write all domain
        j1=jj1(l) ; j2=jj2(l)
        k1=kk1(l) ; k2=kk2(l)
      endif
      
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd = nid*njd

      vartype="Float32"

      open(42,file=file,position="append")

      ! write var

      write(42,'(5A)') '<DataArray type="',trim(vartype),'" Name="',name_field,'"  NumberOfComponents="1" format="ascii">'

      do k=k1,k2-1
          do j=j1,j2-1
            do i=i1,i2-1

            xyz =npc(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

            write(42,'(E16.8)') field(xyz)
          enddo
        enddo
      enddo

      ! write footer

      write(42,'(A)') '</DataArray>'
      close(42)

  end subroutine vtk_writer_r1

  subroutine vtk_writer_r0(vtk_file,field,name_field,l,oi1,oi2,oj1,oj2,ok1,ok2)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      double precision,intent(in) :: field
      integer,intent(in) ::l
      integer,intent(in),optional :: oi1,oi2,oj1,oj2,ok1,ok2
      character(*),intent(in) :: vtk_file,name_field

      integer :: i1,i2,j1,j2,k1,k2
      integer :: nid,njd,nijd,xyz,i,j,k
      character(len=80) :: vartype,file

      if (collect) then
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif
      
      if (present(oi1)) then
        i1=oi1 ; i2=oi2        ! if given, write only a portion of the domain
        j1=oj1 ; j2=oj2
        k1=ok1 ; k2=ok2
      else
        i1=ii1(l) ; i2=ii2(l)  ! if not, write all domain
        j1=jj1(l) ; j2=jj2(l)
        k1=kk1(l) ; k2=kk2(l)
      endif

      vartype="Float32"


      open(42,file=file,position="append")

      ! write var

      write(42,'(5A)') '<DataArray type="',trim(vartype),'" Name="',name_field,'"  NumberOfComponents="1" format="ascii">'

      do k=k1,k2-1
        do j=j1,j2-1
          do i=i1,i2-1
            nid = id2(l)-id1(l)+1
            njd = jd2(l)-jd1(l)+1
            nijd = nid*njd

            xyz =npc(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

            write(42,'(E16.8)') field
          enddo
        enddo
      enddo

      ! write footer

      write(42,'(A)') '</DataArray>'


      close(42)

  end subroutine vtk_writer_r0

  ! THIS IS BUGGY WITH CURRENT GCC VERSION (5.2.0)
  subroutine vtk_writer_poly(vtk_file,x,y,z,field,name_field,l,oi1,oi2,oj1,oj2,ok1,ok2)
      use para_var
      use boundary
      use mod_mpi
      use tools
      implicit none
      double precision,intent(in) :: x(:),y(:),z(:)
      class(*),intent(in) :: field(:)
      integer,intent(in) ::l
      integer,intent(in),optional :: oi1,oi2,oj1,oj2,ok1,ok2
      character(*),intent(in) :: vtk_file,name_field

      integer :: i1,i2,j1,j2,k1,k2
      integer :: nid,njd,nijd,xyz,i,j,k
      character(len=80) :: vartype,file

      if (collect) then
        ncollect=ncollect+1
        call reallocate_s(collection,ncollect)
        collection(ncollect)=vtk_file
        file=trim(collect_dir)//"/"//vtk_file
      else
       file=vtk_file
      endif
      
      if (present(oi1)) then
        i1=oi1 ; i2=oi2        ! if given, write only a portion of the domain
        j1=oj1 ; j2=oj2
        k1=ok1 ; k2=ok2
      else
        i1=ii1(l) ; i2=ii2(l)  ! if not, write all domain
        j1=jj1(l) ; j2=jj2(l)
        k1=kk1(l) ; k2=kk2(l)
      endif

      select type(field)
      type is (integer)
        vartype="Int32"
      type is (real(8))
        vartype="Float32"
      class default
        if(rank==0) write(*,*) 'vtk_writer : Unknown type, must be integer or real'
        stop
      end select

      open(42,file=file,status="replace")

      ! write header

      write(42,'(A)') '<?xml version="1.0"?>'
      write(42,*) '<VTKFile type="StructuredGrid"  version="0.1" byte_order="LittleEndian">'
      write(42,*) '<StructuredGrid WholeExtent="',i1,i2,j1,j2,k1,k2,'">'
      write(42,*) '<Piece Extent="',i1,i2,j1,j2,k1,k2,'">'

      ! write mesh

      write(42,*) "<Points>"
      write(42,*) '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      do k=k1,k2
        do j=j1,j2
          do i=i1,i2

            nid = id2(l)-id1(l)+1
            njd = jd2(l)-jd1(l)+1
            nijd = nid*njd

            xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

            write(42,*) x(xyz),y(xyz),z(xyz)

          enddo
        enddo
      enddo

      write(42,*) '</DataArray>'
      write(42,*) '</Points>'

      ! write var

      write(42,*) '<PointData Scalars="',name_field,'">'
      write(42,*) '<DataArray type="',trim(vartype),'" Name="',name_field,'"  NumberOfComponents="1" format="ascii">'

      if (size(field)==ip21)then ! if it's a field
        do k=k1,k2
          do j=j1,j2
            do i=i1,i2
              nid = id2(l)-id1(l)+1
              njd = jd2(l)-jd1(l)+1
              nijd = nid*njd

              xyz =npn(l)+1+(i -id1(l))+(j -jd1(l))*nid+(k -kd1(l))*nijd

              select type(field)
              type is (integer)
                write(42,*) field(xyz)
              type is (real(8))
                write(42,*) field(xyz)
              class default
                write(*,*) 'vtk_writer : Unknown type, must be integer or real'
                stop
              end select
            enddo
          enddo
        enddo
      elseif  (size(field)==1)then ! if it's a value
        do k=k1,k2
          do j=j1,j2
            do i=i1,i2
              select type(field)
              type is (integer)
                write(42,*) field
              type is (real(8))
                write(42,*) field
              class default
                write(*,*) 'vtk_writer : Unknown type, must be integer or real'
                stop
              end select

            enddo
          enddo
        enddo
      else
        write(*,*) 'vtk_writer : Unknown size of field'
        stop
      endif

      ! write footer

      write(42,*) '</DataArray>'
      write(42,*) "</PointData>"
      write(42,*) '</Piece>'
      write(42,*) '</StructuredGrid>'
      write(42,*) '</VTKFile>'

      close(42)

  end subroutine vtk_writer_poly

end module mod_vtk
