program main
  use iso_fortran_env
  use mpi
  
  implicit none

  ! mpi variables
  integer :: ierr, rank, nprocs

  ! define boundaries of the simulation and the cutoff distance
  double precision, dimension(3) :: lower_boundary = [0,0,0], upper_boundary = [1,1,1]
  double precision :: cutoff = 0.05
  
  ! define positions of the particles
  double precision, dimension(:), allocatable :: globalx, globaly, globalz
  integer :: global_particle_count, local_particle_count

  ! the particles local to the process
  double precision, dimension(:), allocatable :: localx, localy, localz

  ! the particles from other processes
  double precision, dimension(:), allocatable :: foreignx, foreigny, foreignz

  ! stores the pairs
  integer(kind=int64) :: pairs, sum_pairs

  ! mpi boilerplate
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  if (rank == 0) print *, "[TIME] init", elapsed_time()

  ! read the necessary files
  call read_files(globalx, globaly, globalz, global_particle_count, lower_boundary, upper_boundary, cutoff, rank, ierr)
  if (rank == 0) print *, "[TIME] files read", elapsed_time()

  ! calculate the number of particles each process will handle
  if (rank == 0) local_particle_count = global_particle_count / nprocs
  call MPI_BCAST(local_particle_count, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! allocate space in the particle arrays
  allocate(localx(local_particle_count))
  allocate(localy(local_particle_count))
  allocate(localz(local_particle_count))

  allocate(foreignx(local_particle_count))
  allocate(foreigny(local_particle_count))
  allocate(foreignz(local_particle_count))

  ! scatter the particles between the processes
  call MPI_SCATTER(globalx, local_particle_count, MPI_DOUBLE, localx, local_particle_count, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_SCATTER(globaly, local_particle_count, MPI_DOUBLE, localy, local_particle_count, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_SCATTER(globalz, local_particle_count, MPI_DOUBLE, localz, local_particle_count, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

  ! "receive" the local particles into the foreign arrays for comparison
  foreignx = localx
  foreigny = localy
  foreignz = localz

  if (rank == 0) print *, "[TIME] particles distributed", elapsed_time()

  ! do calculations
  pairs = count_pairs(localx, localy, localz, foreignx, foreignz, foreigny, lower_boundary, upper_boundary, cutoff, rank, nprocs, ierr)

  if (rank == 0) then
    print *, "[TIME] counted pairs", elapsed_time()
    print *, "[TIME] total elapsed time", elapsed_time(.true.)
  end if
  ! sum results
  call MPI_REDUCE(pairs, sum_pairs, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! divide by two to counter double counting
  sum_pairs = sum_pairs / 2

  ! print results
  if (rank == 0) print *, " [LOG] total pairs found", sum_pairs

  call MPI_FINALIZE(ierr)

contains

  ! returns the elapsed time since the function was last called
  real(kind=8) function elapsed_time(get_total) result(elapsed)
    implicit none

    ! whether or not to return the running total
    logical, intent(in), optional :: get_total
    
    ! implicit static keeps values between function calls
    real(kind=8) :: start_time = 0.0, end_time = 0.0, running_total = 0.0
    logical :: initialised = .false.
  
    if (present(get_total) .and. initialised) then ! if asking for the total then return it
      elapsed = running_total
    else if (initialised) then ! if the function has already been called once before
      start_time = end_time
      end_time = MPI_WTIME()
      elapsed = end_time - start_time
      running_total = running_total + elapsed
    else ! if this is the first time the function is ran
      start_time = MPI_WTIME()
      initialised = .true.
      elapsed = 0.0
    end if
  end function elapsed_time

  ! reads data from the config and particle data files
  subroutine read_files(posx, posy, posz, &
    particle_count, lower_boundary, upper_boundary, cutoff, &
    rank, ierr)
    implicit none

    ! params
    double precision, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer, intent(inout) :: particle_count

    double precision, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    double precision, intent(inout) :: cutoff

    ! mpi variables
    integer, intent(in) :: rank
    integer, intent(inout) :: ierr
    
    ! temp variables
    integer :: num_particles, seed
    logical :: file_exists

    call read_config(lower_boundary, upper_boundary, cutoff, particle_count, seed, rank, ierr)

    inquire(file="particle_data.txt", exist=file_exists)

    if (file_exists) then
      call read_data(posx, posy, posz, lower_boundary, upper_boundary, rank, ierr)
      particle_count = size(posx)
    else
      call generate_data(posx, posy, posz, seed, particle_count, lower_boundary, upper_boundary, rank, ierr)
    end if
  end subroutine read_files

  ! reads the data from the config file
  subroutine read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed, rank, ierr)
    implicit none

    double precision, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    double precision, intent(inout) :: cutoff
    integer, intent(inout) :: num_particles, seed 

    ! mpi variables
    integer, intent(in) :: rank
    integer, intent(inout) :: ierr

    if (rank == 0) then
      ! open the config file
      open(12, file = "config.txt", status = "old")

      ! read the data from the config file
      read(12, *) num_particles
      read(12, *) seed
      read(12, *) cutoff
      read(12, *) lower_boundary
      read(12, *) upper_boundary

      ! close the config file
      close(12)
    end if

    call MPI_BCAST(num_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cutoff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(lower_boundary, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(upper_boundary, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_config

  ! generates particle data
  subroutine generate_data(posx, posy, posz, seed, num_particles, lower_boundary, upper_boundary, rank, ierr)
    implicit none

    integer, intent(in) :: seed, num_particles
    double precision, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! array to generate the random particles into
    double precision, dimension(:,:), allocatable :: particles 

    ! distance between boundary walls on each axis
    double precision, dimension(3) :: boundary_diff
    
    ! looping varialbes
    integer :: axis, i

    ! mpi variables
    integer, intent(in) :: rank
    integer, intent(inout) :: ierr

    ! random calls require an array of length 8
    integer, dimension(8) :: seed_array
    seed_array = seed

    ! allocate space in all arrays
    allocate(posx(num_particles))
    allocate(posy(num_particles))
    allocate(posz(num_particles))

    if (rank == 0) then
      ! find the distance between the simulation bounds
      do axis = 1, 3
        boundary_diff(axis) = upper_boundary(axis) - lower_boundary(axis)
      end do

      allocate(particles(num_particles,3))

      ! set the seed and then generate the random particles
      call random_seed(put=seed_array)
      call random_number(particles)

      do i = 1, num_particles
        posx(i) = lower_boundary(1) + particles(i, 1) * boundary_diff(1)
        posy(i) = lower_boundary(2) + particles(i, 2) * boundary_diff(2)
        posz(i) = lower_boundary(3) + particles(i, 3) * boundary_diff(3)
      end do

      deallocate(particles)
    end if

    call MPI_BCAST(posx, num_particles, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posy, num_particles, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posz, num_particles, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  end subroutine generate_data
  
  ! reads the file and extracts the particles within the processes domain
  subroutine read_data(posx, posy, posz, &
                       lower_boundary, upper_boundary, &
                       rank, ierr)
    implicit none

    ! variables to be filled with particles and the total number of particles
    double precision, dimension(:), allocatable, intent(inout) :: posx, posy, posz

    ! the bounds of the simulation
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! temporary variables that hold unfiltered data from the file
    integer :: file_length, line

    ! mpi variables
    integer, intent(in) :: rank
    integer, intent(inout) :: ierr

    ! used to iterate over the particles to filter them
    integer :: i

    if (rank == 0) then
      ! define the format of the file and open it
      9 format(f20.17, 4x, f20.17, 4x, f20.17)
      
      open(11, file="particle_data.txt", status="old")
    
      ! get the number of particles in the file
      read(11, *) file_length
    end if

    ! sync
    call MPI_BCAST(file_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! allocate space into the arrays to store the particles
    allocate(posx(file_length))
    allocate(posy(file_length))
    allocate(posz(file_length))

    if (rank == 0) then
      ! read the positions of the particles into the unfiltered array
      do line = 1, file_length
        read(11, 9) posx(line), posy(line), posz(line)
      end do

      ! close file
      close(11)
    end if

    ! sync the values read from the file between the processes
    call MPI_BCAST(posx, file_length, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posy, file_length, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posz, file_length, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_data

  
  ! finds the distance between two points on an axis through a PBC boundary
  double precision function PBC_distance(a, b, lower_boundary, upper_boundary) result(difference)
    implicit none

    ! a and b are two points on an axis within the boundaries
    ! lower and upper boundary are the periodic boundary conditions
    double precision, intent(in) :: a, b, lower_boundary, upper_boundary

    ! the difference will be the sum of the distance between each particle and its closest boundary
    difference = min(upper_boundary - a, a - lower_boundary) + min(upper_boundary - b, b - lower_boundary)
  end function PBC_distance

  ! finds the distance between two points on an axis normally (ignoring PBC)
  double precision function distance(a, b) result(difference)
    implicit none

    ! a and b are two points along an axis
    double precision, intent(in) :: a, b

    ! return the difference between the two points
    difference = abs(a - b)
  end function distance

  ! finds both the PBC distance and the standard distance and returns the smallest value
  double precision function shortest_distance(a, b, lower_boundary, upper_boundary) result(difference)
    implicit none

    ! a and b are two points along an axis
    ! lower and upper boundary are the periodic boundary conditions
    double precision, intent(in) :: a,b, lower_boundary, upper_boundary

    difference = min( &
      distance(a, b), &
      PBC_distance(a, b, lower_boundary, upper_boundary) &
    )
  end function shortest_distance

  ! uses three dimensional pythagoras to find the magnitude of a vector
  double precision function magnitude(vector) result(res)
    implicit none

    ! vector is a array of values in the x, y and z axis
    double precision, dimension(3), intent(in) :: vector

    ! axis stores the axis currently being looped over
    integer :: axis

    ! set the result to 0
    res = 0.0

    ! add the square of each axis' value
    do axis = 1, 3
      res = res + vector(axis) ** 2
    end do

    ! root the result to find the hypotenuse
    res = sqrt(res)
  end function magnitude

  ! finds a given ranks neighbour
  integer function pbc_rank(rank, nprocs) result(pbc)
    implicit none

    integer, intent(in) :: rank, nprocs

    pbc = MOD(rank, nprocs)
    if (pbc < 0) pbc = nprocs + pbc
  end function pbc_rank

  ! checks whether two coordinates are a pair
  logical function check_pair(ax, ay, az, bx, by, bz, lower_boundary, upper_boundary, cutoff) result(is_pair)
    implicit none

    double precision, intent(in) :: ax, ay, az, bx, by, bz, cutoff
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary

    double precision, dimension(3) :: a, b, difference

    integer :: axis

    ! store the values in a format that can be iterated over
    a = [ax, ay, az]
    b = [bx, by, bz]

    do axis = 1, 3
      difference(axis) = shortest_distance( &
        a(axis), b(axis), &
        lower_boundary(axis), upper_boundary(axis))
    end do

    ! check whether the magnitude of the difference is within the cutoff
    is_pair = magnitude(difference) < cutoff
  end function check_pair

  integer(kind=int64) function count_pairs_between_arrays( &
      ax, ay, az, &
      bx, by, bz, &
      lower_boundary, upper_boundary, cutoff, triangular) result(pairs)
    implicit none

    ! the particle arrays to count pairs between
    double precision, dimension(:), allocatable, intent(in) :: ax, ay, az
    double precision, dimension(:), allocatable, intent(in) :: bx, by, bz

    ! the upper and lower boundaries the particles reside in and the max distance between two particles in a pair
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary
    double precision, intent(in) :: cutoff

    ! whether or not to allow j <= i
    logical, intent(in) :: triangular

    ! variables to store the positions of the particles being iterated over and the difference between the two positions
    double precision, dimension(3) :: particle_pos, buffer_pos, difference

    ! variables for iteration
    integer :: i, j_start, j, axis

    pairs = 0

    do i = 1, size(ax)
      do j = 1, size(bx)
        if (triangular .and. j == i) cycle 

        if (check_pair( &
          ax(i), ay(i), az(i), &
          bx(j), by(j), bz(j), &
          lower_boundary, upper_boundary, cutoff &
        )) then
          pairs = pairs + 1
        end if
      end do
    end do
  end function count_pairs_between_arrays

  integer(kind=int64) function count_pairs(localx, localy, localz, foreignx, foreigny, foreignz, lower_boundary, upper_boundary, cutoff, rank, nprocs, ierr) result(pairs)
    implicit none

    ! the original particles local to the process
    double precision, dimension(:), allocatable, intent(in) :: localx, localy, localz

    ! the particles to compare against (foreign to the process)
    double precision, dimension(:), allocatable, intent(in) :: foreignx, foreigny, foreignz

    ! the upper and lower boundaries the particles reside in and the max distance between two particles in a pair
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary
    double precision, intent(in) :: cutoff

    ! arrays to store values to be sent to the next process
    double precision, dimension(:), allocatable :: sendx, sendy, sendz

    ! mpi variables
    integer, intent(in) :: rank, nprocs
    integer, intent(inout) :: ierr
    integer :: status0(MPI_STATUS_SIZE)
    integer :: next_rank, prev_rank

    ! variables used for iteration
    integer :: i

    ! allocate space in the send buffers
    allocate(sendx(size(foreignx)))
    allocate(sendy(size(foreignx)))
    allocate(sendz(size(foreignx)))

    ! start pairs at 0
    pairs = 0

    ! find the neighbour to the right and left of the current rank
    next_rank = pbc_rank(rank + 1, nprocs)
    prev_rank = pbc_rank(rank - 1, nprocs)

    do i = 1, nprocs
      ! find pairs between the original particles and the buffer
      pairs = pairs + count_pairs_between_arrays( &
        localx, localy, localz, &
        foreignx, foreigny, foreignz, &
        lower_boundary, upper_boundary, cutoff, &
        i == 1)

      ! prepare the buffer to be sent to the next process
      sendx = foreignx
      sendy = foreigny
      sendz = foreignz

      ! send the data to the next process and rewrite the buffer with the received data
      call MPI_SENDRECV(sendx, size(sendx), MPI_DOUBLE, next_rank, rank, &
                        foreignx, size(sendx), MPI_DOUBLE, prev_rank, prev_rank, &
                        MPI_COMM_WORLD, status0, ierr)

      call MPI_SENDRECV(sendy, size(sendy), MPI_DOUBLE, next_rank, rank, &
                        foreigny, size(sendy), MPI_DOUBLE, prev_rank, prev_rank, &
                        MPI_COMM_WORLD, status0, ierr)

      call MPI_SENDRECV(sendz, size(sendz), MPI_DOUBLE, next_rank, rank, &
                        foreignz, size(sendz), MPI_DOUBLE, prev_rank, prev_rank, &
                        MPI_COMM_WORLD, status0, ierr)
    end do

    deallocate(sendx)
    deallocate(sendy)
    deallocate(sendz)
  end function count_pairs
end program main
