    !----------------------------------------------------------------
    ! Convert ramses binary fortran output to hdf5 type outputs.
    ! For now, it only converts dark matter particles into a simple
    ! hdf5 file, but it can be extended to be used as a converter 
    ! for ramses output to other simulation type outputs, like
    ! gadget or swift. It could also be extended to work in parallel.
    ! If you need anything of this sort, contact me via 
    ! mladen.ivkovic@hotmail.com
    !
    ! Usage:
    !     ramses2hdf5 [options] path/to/output_XXXXX
    !
    ! options for [options] :
    !     -d | --dmo )
    !        work for dark matter particles only. Assumes all particles
    !        have identical mass. If you have variable particle mass,
    !        you need to use the -vm or --variable-mass flag as well.
    !        Is the default setting if you specify no options.
    !
    !     -g | --gadget )
    !        create a GADGET compatible hdf5 type file.
    !        NOT IMPLEMENTED YET
    !
    !     --galaxies )
    !        also write galaxies and orphan galaxies as created by
    !        the mergertree patch as individual groups in the hdf5 file.
    !
    !     -sw | --swift )
    !        create a SWIFT compatible hdf5 type file.
    !        NOT IMPLEMENTED YET
    !
    !     -vm | --variable-mass )
    !        assume DM particles have variable masses, and create a
    !        group for each particle mass.
    !        NOT IMPLEMENTED YET
    !
    ! The program will create an output .h5 file in the output_XXXXX
    ! directory.
    !----------------------------------------------------------------


program ramses2hdf5

  use hdf5

  implicit none

  integer, parameter :: dp = kind(1.d0)

  ! for which output case to work
  logical :: dmo      = .true.   ! use DM particles only
  logical :: gadget   = .false.   ! create gadget compatible file
  logical :: galaxies = .true.   ! use galaxies and orphans from mergertree patch
  logical :: swift    = .false.   ! create swift compatible file
  logical :: varmass  = .false.   ! assume variable DM particle masses

  character(len=100) :: sourcedir = ""  ! output directory 
  character(len=100) :: infofile = ""   ! info_XXXXX.txt file
  character(len=100) :: outputfile = "" ! info_XXXXX.txt file

  ! snapshot data
  integer :: ncpu
  integer :: nstep_coarse
  integer :: snapshot

  real(dp) :: boxlen
  real(dp) :: t
  real(dp) :: aexp
  real(dp) :: H0
  real(dp) :: omega_m
  real(dp) :: omega_l
  real(dp) :: omega_k
  real(dp) :: omega_b
  real(dp) :: unit_l
  real(dp) :: unit_d
  real(dp) :: unit_t

  ! particle counters
  integer :: nparttot = 0
  integer :: ngalaxies = 0
  integer :: norphans = 0

  ! global hdf5 stuff
  integer (HID_T) :: file_id

  ! misc
  integer :: error


  call get_cmdlineargs()
  call read_infofile()

  ! Fire up hdf5 and create file
  call h5open_f(error)
  call h5fcreate_f(outputfile, H5F_ACC_TRUNC_F, file_id, error)


  ! if (dmo) call write_dmo
  call write_dmo
  call write_galaxies
  call write_metadata()

  ! close up hdf5 and file
  call h5fclose_f(file_id, error)
  call h5close_f(error)

  write(*, '(A)') " Finished."


contains

  !==========================================
  subroutine write_galaxies()
  !==========================================
  ! Read in and write galaxy particles
  ! to the hdf5 file.
  ! Ignores particle masses, just writes 
  ! positions.
  !------------------------------------------

    implicit none

    character(len=13),  parameter :: galgroupname = "Galaxies"        ! Group name
    character(len=13),  parameter :: orphgroupname = "Orphans"        ! Group name
    character(len=13),  parameter :: datasetname = "positions"     ! dataset name

    integer(HID_T) :: galgroup_id       ! Group identifier
    integer(HID_T) :: orphgroup_id      ! Group identifier
    integer(HID_T) :: galdataset_id     ! Dataset identifier
    integer(HID_T) :: orphdataset_id    ! Dataset identifier
    integer(HID_T) :: galdataspace_id   ! Data space identifier
    integer(HID_T) :: orphdataspace_id  ! Data space identifier
    integer(HID_T) :: galmemspace_id    ! memory space identifier
    integer(HID_T) :: orphmemspace_id   ! memory space identifier
    integer(HID_T) :: crp_list      ! dataset creation property identifier

    ! Dataset dimensions in the file
    integer(HSIZE_T), dimension(1:2) :: dset_dims ! dataset dimensions at creation time
    integer(HSIZE_T), dimension(1:2) :: chunk_dims
    integer(HSIZE_T), dimension(1:2) :: memory_dims ! dimensions of memory space; In our case = count(:, :)

    ! maximum dimensions
    integer(HSIZE_T), dimension(1:2) :: maxdims
    integer(HSIZE_T), dimension(1:2) :: galoffset, galcount
    integer(HSIZE_T), dimension(1:2) :: orphoffset, orphcount

    ! data
    real(dp), dimension(:, :), allocatable :: xgals, xorph

    integer(HSIZE_T), dimension(1:2) :: data_dims
    integer(HSIZE_T), dimension(1:2) :: galsize, orphsize

    integer(HSIZE_T) :: i, j, r
    integer(HSIZE_T) :: rows, columns, ngals, norph
    integer(HSIZE_T) :: zero = 0
    integer :: error, rank=2

    error = 0 


    ! Create a group in the file.
    call h5gcreate_f(file_id, galgroupname, galgroup_id, error)
    call h5gcreate_f(file_id, orphgroupname, orphgroup_id, error)

    ! Create the data space with unlimited dimensions.
    ! create some junk data to initialize dataset
    rows = 3
    columns = 1
    dset_dims = (/rows, columns/)
    maxdims = (/rows, H5S_UNLIMITED_F/)
    call h5screate_simple_f(rank, dset_dims, galdataspace_id, error, maxdims)
    call h5screate_simple_f(rank, dset_dims, orphdataspace_id, error, maxdims)

    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    chunk_dims = (/rows, columns/)
    call h5pset_chunk_f(crp_list, rank, chunk_dims, error)

    ! Create a dataset in group with default properties.
    CALL h5dcreate_f(galgroup_id, datasetname, H5T_NATIVE_DOUBLE, galdataspace_id, &
         galdataset_id, error, crp_list)
    CALL h5dcreate_f(orphgroup_id, datasetname, H5T_NATIVE_DOUBLE, orphdataspace_id, &
         orphdataset_id, error, crp_list)

    call h5sclose_f(galdataspace_id, error)
    call h5sclose_f(orphdataspace_id, error)

    ! fill initial data array (with junk) and write it to dataset
    data_dims = (/rows, columns/)
    allocate(xorph(1:rows, 1:columns))
    do i = 1, data_dims(1)
      do j = 1, data_dims(2)
        xorph(i,j) = float(-1)
      enddo
    enddo

    call h5dwrite_f(galdataset_id, H5T_NATIVE_DOUBLE, xorph, data_dims, error)
    call h5dwrite_f(orphdataset_id, H5T_NATIVE_DOUBLE, xorph, data_dims, error)
    deallocate(xorph)

    galoffset(1:2) = (/zero, zero/) ! set initial offset to zero; overwrite initial junk
    galcount(1:2) = (/zero, zero/)
    galsize(1:2) = (/rows, zero/)
    orphoffset(1:2) = (/zero, zero/)
    orphcount(1:2) = (/zero, zero/)
    orphsize(1:2) = (/rows, zero/)

    do r = 1, ncpu

      call read_galaxy_positions(r, xgals, xorph, ngals, norph)

      ! FIRST DO GALAXIES
      data_dims = (/rows, ngals/)

      ! extend the dataset
      galsize(2) = galsize(2) + ngals
      call h5dset_extent_f(galdataset_id, galsize, error)

      ! Create the memory space for the selection
      galcount = (/rows, ngals/)
      memory_dims = galcount
      call h5screate_simple_f(rank, memory_dims, galmemspace_id, error)

      ! create data and write to extended part of dataset
      call h5dget_space_f(galdataset_id, galdataspace_id, error)
      call h5sselect_hyperslab_f(galdataspace_id, H5S_SELECT_SET_F, galoffset, galcount, error)
      ! update offset for next round
      galoffset(2) = galoffset(2) + ngals

      call h5dwrite_f(galdataset_id, H5T_NATIVE_DOUBLE, xgals, data_dims, error, galmemspace_id, galdataspace_id)
      deallocate(xgals)


      ! NOW DO ORPHANS
      data_dims = (/rows, norph/)

      ! extend the dataset
      orphsize(2) = orphsize(2) + norph
      call h5dset_extent_f(orphdataset_id, orphsize, error)

      ! Create the memory space for the selection
      orphcount = (/rows, norph/)
      memory_dims = orphcount
      call h5screate_simple_f(rank, memory_dims, orphmemspace_id, error)

      ! create data and write to extended part of dataset
      call h5dget_space_f(orphdataset_id, orphdataspace_id, error)
      call h5sselect_hyperslab_f(orphdataspace_id, H5S_SELECT_SET_F, orphoffset, orphcount, error)
      ! update offset for next round
      orphoffset(2) = orphoffset(2) + norph

      call h5dwrite_f(orphdataset_id, H5T_NATIVE_DOUBLE, xorph, data_dims, error, orphmemspace_id, orphdataspace_id)
      deallocate(xorph)
    enddo

    ! Close dataspace, dataset, and group
    call h5dclose_f(galdataset_id, error)
    call h5sclose_f(galmemspace_id, error)
    call h5gclose_f(galgroup_id, error)
    call h5dclose_f(orphdataset_id, error)
    call h5sclose_f(orphmemspace_id, error)
    call h5gclose_f(orphgroup_id, error)


  end subroutine write_galaxies





  !==========================================
  subroutine write_dmo()
  !==========================================
  ! Read in and write dark matter particles
  ! to the hdf5 file.
  ! Ignores particle masses, just writes 
  ! positions.
  !------------------------------------------

    implicit none

    character(len=13),  parameter :: groupname = "DMparticles"   ! Group name
    character(len=13),  parameter :: datasetname = "positions"     ! dataset name

    integer(HID_T) :: group_id      ! Group identifier
    integer(HID_T) :: dataset_id    ! Dataset identifier
    integer(HID_T) :: dataspace_id  ! Data space identifier
    integer(HID_T) :: memspace_id   ! memory space identifier
    integer(HID_T) :: crp_list      ! dataset creation property identifier

    ! Dataset dimensions in the file
    integer(HSIZE_T), dimension(1:2) :: dset_dims ! dataset dimensions at creation time
    integer(HSIZE_T), dimension(1:2) :: chunk_dims
    integer(HSIZE_T), dimension(1:2) :: memory_dims ! dimensions of memory space; In our case = count(:, :)

    ! maximum dimensions
    integer(HSIZE_T), dimension(1:2) :: maxdims
    integer(HSIZE_T), dimension(1:2) :: offset, count

    ! data
    real(dp), dimension(:, :), allocatable :: x

    integer(HSIZE_T), dimension(1:2) :: data_dims
    integer(HSIZE_T), dimension(1:2) :: size

    integer(HSIZE_T) :: i, j, r
    integer(HSIZE_T) :: rows, columns, nparts
    integer(HSIZE_T) :: zero = 0
    integer :: error, rank=2

    error = 0 


    ! Create a group in the file.
    call h5gcreate_f(file_id, groupname, group_id, error)

    ! Create the data space with unlimited dimensions.
    ! create some junk data to initialize dataset
    rows = 3
    columns = 1
    dset_dims = (/rows, columns/)
    maxdims = (/rows, H5S_UNLIMITED_F/)
    call h5screate_simple_f(rank, dset_dims, dataspace_id, error, maxdims)

    ! Modify dataset creation properties, i.e. enable chunking
    call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    chunk_dims = (/rows, columns/)
    call h5pset_chunk_f(crp_list, rank, chunk_dims, error)

    ! Create a dataset in group with default properties.
    CALL h5dcreate_f(group_id, datasetname, H5T_NATIVE_DOUBLE, dataspace_id, &
         dataset_id, error, crp_list)

    call h5sclose_f(dataspace_id, error)

    ! fill initial data array (with junk) and write it to dataset
    data_dims = (/rows, columns/)
    allocate(x(1:rows, 1:columns))
    do i = 1, data_dims(1)
      do j = 1, data_dims(2)
        x(i,j) = float(-1)
      enddo
    enddo

    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x, data_dims, error)
    deallocate(x)

    offset(1:2) = (/zero, zero/) ! set initial offset to zero; overwrite initial junk
    count(1:2) = (/zero, zero/)
    size(1:2) = (/rows, zero/)

    do r = 1, ncpu

      call read_particle_positions(r, x, nparts)
      data_dims = (/rows, nparts/)

      ! extend the dataset
      size(2) = size(2) + nparts
      call h5dset_extent_f(dataset_id, size, error)

      ! Create the memory space for the selection
      count = (/rows, nparts/)
      memory_dims = count
      call h5screate_simple_f(rank, memory_dims, memspace_id, error)

      ! create data and write to extended part of dataset
      call h5dget_space_f(dataset_id, dataspace_id, error)
      call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
      ! update offset for next round
      offset(2) = offset(2) + nparts

      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, x, data_dims, error, memspace_id, dataspace_id)
      deallocate(x)

    enddo

    ! Close dataspace, dataset, and group
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(memspace_id, error)
    call h5gclose_f(group_id, error)


  end subroutine write_dmo






  !============================================================
  subroutine read_particle_positions(cpu, x, nparts) 
  !============================================================
  ! Reads in particle positions from part_XXXXX.outYYYYY files
  !------------------------------------------------------------

    implicit none
    integer(HSIZE_T), intent(in)           :: cpu
    real(dp), allocatable, dimension(:, :) :: x
    integer(HSIZE_T), intent(out)          :: nparts

    character(len=200) :: fname
    character(len=5)   :: ids, sns

    integer :: junk, nparts_temp

    logical :: exists

    call title(int(cpu), ids)
    call title(snapshot, sns)

    fname = TRIM(sourcedir)//'/part_'//TRIM(sns)//'.out'//TRIM(ids)
    write(*,*) "Working on file ", fname 
    inquire(file=fname, exist=exists)
    if (.not.exists) then
      write(*,*) "Can't find file ", fname
      write(*,*) "exiting"
      stop
    endif

    open(666, file=fname, form='unformatted')
    read(666) junk
    read(666) junk
    read(666) nparts_temp

    read(666) junk
    read(666) junk
    read(666) junk
    read(666) junk
    read(666) junk

    allocate(x(1:3, 1:nparts_temp))
    read(666) x(1, :)
    read(666) x(2, :)
    read(666) x(3, :)

    close(666)

    nparttot = nparttot + nparts_temp
    nparts = int(nparts_temp, kind=HSIZE_T)


  end subroutine read_particle_positions






  !=================================================================
  subroutine read_galaxy_positions(cpu, xgal, xorph, ngals, norph) 
  !=================================================================
  ! Reads in galaxy positions from part_XXXXX.outYYYYY files
  !-----------------------------------------------------------------

    implicit none
    integer(HSIZE_T), intent(in)           :: cpu
    real(dp), allocatable, dimension(:, :) :: xorph
    real(dp), allocatable, dimension(:, :) :: xgal
    integer(HSIZE_T), intent(out)          :: ngals
    integer(HSIZE_T), intent(out)          :: norph

    real(dp), allocatable, dimension(:, :) :: xall
    integer, allocatable, dimension(:) :: clmpid

    character(len=200) :: fname
    character(len=5)   :: ids, sns

    logical :: file_exists    
    integer :: ngals_tot, ngals_temp, norph_temp, io, i, j, igal, iorph
    real(dp):: junk_real

    call title(int(cpu), ids)
    call title(snapshot, sns)

    fname = TRIM(sourcedir)//'/galaxies_'//TRIM(sns)//'.txt'//TRIM(ids)
    write(*,*) "Working on file ", fname 

    ngals_tot = 0

    inquire(file=fname, exist=file_exists)
    if (file_exists) then
      open(unit=600, file=fname)
      read(600, *) ! skip header
      do 
        read(600,*, iostat=io)
        if (io/=0) exit
        ngals_tot = ngals_tot + 1
      enddo
      close(600)
    else
      write(*,*) "Didn't find file ", fname
    endif


    ngals_temp = 0
    norph_temp = 0
  
    allocate(xall(1:3, 1:ngals_tot))
    allocate(clmpid(1:ngals_tot))

    open(unit=600, file=fname)
    read(600,*) ! skip header
    do i=1, ngals_tot
      read(600,'(I20,x,4(E20.12,x))') clmpid(i), junk_real, xall(1, i), xall(2, i), xall(3, i) 
      if (clmpid(i)>0) then
        ngals_temp = ngals_temp + 1
      else
        norph_temp = norph_temp + 1
      endif
    enddo
    close(600)


    allocate(xgal(1:3, 1:ngals_temp))
    allocate(xorph(1:3, 1:norph_temp))

    igal = 1
    iorph = 1
    do i = 1, ngals_tot
      if (clmpid(i)>0) then
        do j=1, 3
          xgal(j,igal) = xall(j, i)
        enddo
        igal = igal + 1
      else
        do j=1, 3
          xorph(j,iorph) = xall(j, i)
        enddo
        iorph = iorph + 1
      endif
    enddo



    ngalaxies = ngalaxies + ngals_temp
    ngals = int(ngals_temp, kind=HSIZE_T)
    norphans = norphans + norph_temp
    norph = int(norph_temp, kind=HSIZE_T)

    deallocate(xall, clmpid)



  end subroutine read_galaxy_positions




  !=========================================
  subroutine write_metadata()
  !=========================================
    ! write snapshot metadata as attributes
    ! in a group in the file
    !---------------------------------------

    implicit none
    character(len=8),  parameter :: groupname = "Header"     ! Group name; needs to be Header for VR
    integer(HID_T)               :: group_id                 ! Group identifier
    integer                      :: error                    ! Error flag

    ! Create a group in the file.
    call h5gcreate_f(file_id, groupname, group_id, error)

    call add_group_attribute_int(group_id, "snapshot", snapshot)
    call add_group_attribute_int(group_id, "nstep_coarse", nstep_coarse)
    call add_group_attribute_int(group_id, "ncpu", ncpu)
    call add_group_attribute_int(group_id, "nparttot", nparttot)
    if (galaxies) then
      call add_group_attribute_int(group_id, "ngalaxies", ngalaxies)
      call add_group_attribute_int(group_id, "norphans", norphans)
    endif

    call add_group_attribute_real(group_id, "BoxSize", boxlen) ! Needs to be BoxSize for VR
    call add_group_attribute_real(group_id, "t", t)
    call add_group_attribute_real(group_id, "aexp", aexp)
    call add_group_attribute_real(group_id, "H0", H0)
    call add_group_attribute_real(group_id, "omega_m", omega_m)
    call add_group_attribute_real(group_id, "omega_l", omega_l)
    call add_group_attribute_real(group_id, "omega_b", omega_b)
    call add_group_attribute_real(group_id, "omega_k", omega_k)
    call add_group_attribute_real(group_id, "unit_l", unit_l)
    call add_group_attribute_real(group_id, "unit_d", unit_d)
    call add_group_attribute_real(group_id, "unit_t", unit_t)

    ! Close the group.
    CALL h5gclose_f(group_id, error)


  end subroutine write_metadata





  !=============================================================
  subroutine add_group_attribute_int(group_id, aname, aval)
  !=============================================================
    ! add a integer group attribute
    !-----------------------------------------------------------

    implicit none

    integer(HID_T), intent(in)   :: group_id
    character(len=*), intent(in) :: aname
    integer, intent(in)          :: aval

    integer(HID_T) :: attr_id       ! Attribute identifier
    integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier

    integer(HSIZE_T), dimension(1)  :: adims = (/1/)  ! Attribute dimension
    integer                         :: arank = 1      ! Attribure rank

    CALL h5screate_simple_f(arank, adims, aspace_id, error)
    CALL h5acreate_f(group_id, aname, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
    CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, aval, adims, error)
    CALL h5aclose_f(attr_id, error)
    CALL h5sclose_f(aspace_id, error)

  end subroutine add_group_attribute_int




  !=============================================================
  subroutine add_group_attribute_real(group_id, aname, aval)
  !=============================================================
    ! add a float group attribute
    !-----------------------------------------------------------

    implicit none

    integer(HID_T), intent(in)   :: group_id
    character(len=*), intent(in) :: aname
    real(dp), intent(in)         :: aval

    integer(HID_T) :: attr_id       ! Attribute identifier
    integer(HID_T) :: aspace_id     ! Attribute Dataspace identifier

    integer(HSIZE_T), dimension(1)  :: adims = (/1/)  ! Attribute dimension
    integer                         :: arank = 1      ! Attribure rank

    CALL h5screate_simple_f(arank, adims, aspace_id, error)
    CALL h5acreate_f(group_id, aname, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, aval, adims, error)
    CALL h5aclose_f(attr_id, error)
    CALL h5sclose_f(aspace_id, error)

  end subroutine add_group_attribute_real




  !================================ 
  subroutine read_infofile()
  !================================ 

    implicit none
    character(len=13) :: junk_string
    integer  :: junk_int

    open(unit=666, file=infofile, form='formatted')
    read(666, '(A13,I11)') junk_string, ncpu
    read(666, '(A13,I11)') junk_string, junk_int
    read(666, '(A13,I11)') junk_string, junk_int 
    read(666, '(A13,I11)') junk_string, junk_int 
    read(666, '(A13,I11)') junk_string, junk_int 
    read(666, '(A13,I11)') junk_string, nstep_coarse
    read(666, *)
    read(666, '(A13,E23.15)') junk_string, boxlen
    read(666, '(A13,E23.15)') junk_string, t
    read(666, '(A13,E23.15)') junk_string, aexp
    read(666, '(A13,E23.15)') junk_string, H0
    read(666, '(A13,E23.15)') junk_string, omega_m
    read(666, '(A13,E23.15)') junk_string, omega_l
    read(666, '(A13,E23.15)') junk_string, omega_k
    read(666, '(A13,E23.15)') junk_string, omega_b
    read(666, '(A13,E23.15)') junk_string, unit_l
    read(666, '(A13,E23.15)') junk_string, unit_d
    read(666, '(A13,E23.15)') junk_string, unit_t
    close(666)

  end subroutine read_infofile



  !======================================
  subroutine get_cmdlineargs()
  !======================================
    ! read and sort out cmdlineargs.
    ! also generate and check sourcedir
    ! and infofile.
    !------------------------------------
  
    implicit none

    integer            :: i, strlen
    character(len=100) :: arg
    character(len=15)  :: ftype
    logical            :: dir_ex, infofile_ex

    if (command_argument_count()==0) then
      call print_help()
    else
      do i=1, command_argument_count()
        call get_command_argument(i, arg)
        select case(trim(arg))
          case ('-d', '--dmo')
            dmo = .true.
          case ('-g', '--gadget')
            gadget=.true.
          case ('--galaxies')
            galaxies = .true.
          case ('-sw', '--swift')
            swift = .true.
          case ('-vm', '--variable-massses')
            varmass = .true.
          case default
            sourcedir = trim(arg)
        end select 
      enddo
    endif



    ! Now do some preliminary checks

    if (gadget) then
      write(*,*) "Conversion to gadget format not implemented yet."
      stop
    endif
    if (swift) then
      write(*,*) "Conversion to swift format not implemented yet."
      stop
    endif
    if (varmass) then
      write(*,*) "Variable masses not implemented yet."
      stop
    endif

    ! if nothing else is specified, work for DM only
    if ((.not.gadget) .and. (.not.swift) .and. (.not.galaxies)) dmo=.true.

    ! check whether output_XXXXX directory is given
    strlen = len(trim(sourcedir))
    if (strlen==0) then
      write(*,*) "ERROR: no output directory given. Exiting."
      call print_help()
    else
      inquire(file=TRIM(sourcedir), exist=dir_ex)
      if (.not.dir_ex) then
        write(*,*) "Couldn't find directory", TRIM(sourcedir), "; Exiting."
        stop
      else

        ! get ftype to also generate output file type at this point
        if (gadget)           ftype='gadget'
        if (swift)            ftype='swift'
        if (dmo)              ftype='dmo'
        if (galaxies)         ftype='galaxies'
        if (dmo.and.galaxies) ftype='dmo_galaxies'

        if(sourcedir(strlen:strlen)=='/') then
          infofile = TRIM(sourcedir)//'info_'//sourcedir(strlen-5:strlen-1)//'.txt'
          outputfile = TRIM(sourcedir)//'output-'//TRIM(ftype)//'_'//sourcedir(strlen-5:strlen-1)//'.hdf5'
          read(sourcedir(strlen-5:strlen-1), '(I5)') snapshot
        else
          infofile = TRIM(sourcedir)//'/info_'//sourcedir(strlen-4:strlen)//'.txt'
          outputfile = TRIM(sourcedir)//'/output-'//TRIM(ftype)//'_'//sourcedir(strlen-4:strlen)//'.hdf5'
          read(sourcedir(strlen-4:strlen), '(I5)') snapshot
        endif
        inquire(file=TRIM(infofile), exist=infofile_ex)
        if (.not. infofile_ex) then
          write(*,*) "Couldn't find info file ", infofile
          write(*,*) "Exiting"
          stop
        endif
      endif
    endif

    write(*,*) "Results will be written to ", outputfile


  end subroutine get_cmdlineargs



  !======================================
  subroutine print_help()
  !======================================
    ! print the help message and quit.
    !------------------------------------

    implicit none

    write(*,*) "Help message will be done soon"
    stop

  end subroutine print_help

  !======================================
  subroutine title(n,nchar)
  !======================================
    integer::n
    character(5)::nchar

    character(1)::nchar1
    character(2)::nchar2
    character(3)::nchar3
    character(4)::nchar4
    character(5)::nchar5

    if(n.ge.10000)then
       write(nchar5,'(i5)') n
       nchar = nchar5
    elseif(n.ge.1000)then
       write(nchar4,'(i4)') n
       nchar = '0'//nchar4
    elseif(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = '00'//nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '000'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '0000'//nchar1
    endif

  end subroutine title




end program ramses2hdf5
