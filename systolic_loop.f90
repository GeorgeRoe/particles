program main
  use iso_fortran_env
  use mpi
  
  implicit none

  ! mpi variables
  integer :: ierr, rank, nprocs

  ! define boundaries of the simulation and the cutoff distance
  real, dimension(3) :: lower_boundary = [0,0,0], upper_boundary = [1,1,1]
  real :: cutoff = 0.5
  
  ! temporary variables in which the result of the file read is stored
  real, dimension(:), allocatable :: readx, ready, readz
  integer :: file_length

  ! define positions of the particles
  real, dimension(:), allocatable :: posx, posy, posz
  integer :: particle_count

  ! the buffers in which the particles are stored
  real, dimension(:), allocatable :: bufferx, buffery, bufferz

  ! stores the pairs
  integer(kind=int64) :: pairs, sum_pairs

  ! mpi boilerplate
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  ! read the file
  if (rank == 0) call read_file(readx, readz, ready, file_length)

  ! calculate the number of particles each process will handle
  if (rank == 0) particle_count = file_length / nprocs
  call MPI_BCAST(particle_count, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! allocate space in the particle arrays
  allocate(posx(particle_count))
  allocate(posy(particle_count))
  allocate(posz(particle_count))

  allocate(bufferx(particle_count))
  allocate(buffery(particle_count))
  allocate(bufferz(particle_count))

  ! scatter the particles between the processes
  call MPI_SCATTER(readx, particle_count, MPI_REAL, posx, particle_count, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_SCATTER(ready, particle_count, MPI_REAL, posy, particle_count, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_SCATTER(readz, particle_count, MPI_REAL, posz, particle_count, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

  ! deallocate the space for reading the file
  if (rank == 0) then
    deallocate(readx)
    deallocate(ready)
    deallocate(readz)
  end if

  ! store the current particles in the buffer to be compared
  bufferx = posx
  buffery = posy
  bufferz = posz

  ! do calculations
  pairs = count_pairs(posx, posy, posz, bufferx, bufferz, buffery, lower_boundary, upper_boundary, cutoff, rank, nprocs, ierr)
  ! sum results
  call MPI_REDUCE(pairs, sum_pairs, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! divide by two to counter double counting
  sum_pairs = sum_pairs / 2

  ! print results
  if (rank == 0) print *, "total pairs found", sum_pairs

  call MPI_FINALIZE(ierr)

contains

  ! reads particle data from the particle data file
  subroutine read_file(posx, posy, posz, file_length)
    implicit none

    ! files are formatted as such:
    ! the first line contains how many rows worth of data the file stores
    ! the subsequent lines statustore data in columns for the x, y and z positions of the particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer, intent(inout) :: file_length

    ! open the file
    open(11, file = "particle_data.txt", status = "old")

    ! read how many lines of data the file holds and allocate the space
    read(11, *) file_length

    allocate(posx(file_length))
    allocate(posy(file_length))
    allocate(posz(file_length))

    ! define the format of the columnar data and read it into the variables
    9 format(f20.17, 4x, f20.17, 4x, f20.17)
    read(11, 9) posx, posy, posz

    close(11)
  end subroutine read_file
  
  ! finds the distance between two points on an axis through a PBC boundary
  real function PBC_distance(a, b, lower_boundary, upper_boundary) result(difference)
    implicit none

    ! a and b are two points on an axis within the boundaries
    ! lower and upper boundary are the periodic boundary conditions
    real, intent(in) :: a, b, lower_boundary, upper_boundary

    ! the difference will be the sum of the distance between each particle and its closest boundary
    difference = min(upper_boundary - a, a - lower_boundary) + min(upper_boundary - b, b - lower_boundary)
  end function PBC_distance

  ! finds the distance between two points on an axis normally (ignoring PBC)
  real function distance(a, b) result(difference)
    implicit none

    ! a and b are two points along an axis
    real, intent(in) :: a, b

    ! return the difference between the two points
    difference = abs(a - b)
  end function distance

  ! finds both the PBC distance and the standard distance and returns the smallest value
  real function shortest_distance(a, b, lower_boundary, upper_boundary) result(difference)
    implicit none

    ! a and b are two points along an axis
    ! lower and upper boundary are the periodic boundary conditions
    real, intent(in) :: a,b, lower_boundary, upper_boundary

    difference = min( &
      distance(a, b), &
      PBC_distance(a, b, lower_boundary, upper_boundary) &
    )
  end function shortest_distance

  ! uses three dimensional pythagoras to find the magnitude of a vector
  real function magnitude(vector) result(res)
    implicit none

    ! vector is a array of values in the x, y and z axis
    real, dimension(3), intent(in) :: vector

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

    real, intent(in) :: ax, ay, az, bx, by, bz, cutoff
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    real, dimension(3) :: a, b, difference

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

  integer(kind=int64) function count_pairs_against_buffer( &
      posx, posy, posz, &
      bufferx, buffery, bufferz, &
      lower_boundary, upper_boundary, cutoff, triangular) result(pairs)
    implicit none

    ! the original particles recieved by the processor
    real, dimension(:), allocatable, intent(in) :: posx, posy, posz

    ! the particles currently stored in the buffer
    real, dimension(:), allocatable, intent(in) :: bufferx, buffery, bufferz

    ! the upper and lower boundaries the particles reside in and the max distance between two particles in a pair
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary
    real, intent(in) :: cutoff

    ! whether or not to allow j <= i
    logical, intent(in) :: triangular

    ! variables to store the positions of the particles being iterated over and the difference between the two positions
    real, dimension(3) :: particle_pos, buffer_pos, difference

    ! variables for iteration
    integer :: i, j_start, j, axis

    pairs = 0

    do i = 1, size(posx)
      do j = 1, size(bufferx)
        if (triangular .and. j == i) cycle 

        if (check_pair( &
          posx(i), posy(i), posz(i), &
          bufferx(j), buffery(j), bufferz(j), &
          lower_boundary, upper_boundary, cutoff &
        )) then
          pairs = pairs + 1
        end if
      end do
    end do
  end function count_pairs_against_buffer

  integer(kind=int64) function count_pairs(posx, posy, posz, bufferx, bufferz, buffery, lower_boundary, upper_boundary, cutoff, rank, nprocs, ierr) result(pairs)
    implicit none

    ! the original particles recieved by the processor
    real, dimension(:), allocatable, intent(in) :: posx, posy, posz

    ! the particles currently stored in the buffer
    real, dimension(:), allocatable, intent(in) :: bufferx, buffery, bufferz

    ! the upper and lower boundaries the particles reside in and the max distance between two particles in a pair
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary
    real, intent(in) :: cutoff

    ! arrays to store values to be sent to the next process
    real, dimension(:), allocatable :: sendx, sendy, sendz

    ! mpi variables
    integer, intent(in) :: rank, nprocs
    integer, intent(inout) :: ierr
    integer :: status0(MPI_STATUS_SIZE)
    integer :: next_rank, prev_rank

    ! variables used for iteration
    integer :: i
    integer(kind=int64) :: current_pairs

    ! allocate space in the send buffers
    allocate(sendx(size(bufferx)))
    allocate(sendy(size(bufferx)))
    allocate(sendz(size(bufferx)))

    ! start pairs at 0
    pairs = 0

    ! find the neighbour to the right and left of the current rank
    next_rank = pbc_rank(rank + 1, nprocs)
    prev_rank = pbc_rank(rank - 1, nprocs)

    do i = 1, nprocs
      ! find pairs between the original particles and the buffer
      current_pairs = count_pairs_against_buffer( &
        posx, posy, posz, &
        bufferx, buffery, bufferz, &
        lower_boundary, upper_boundary, cutoff, &
        i == 1)
      pairs = pairs + current_pairs

      ! prepare the buffer to be sent to the next process
      sendx = bufferx
      sendy = buffery
      sendz = bufferz

      ! send the data to the next process and rewrite the buffer with the received data
      call MPI_SEND(sendx, size(bufferx), MPI_REAL, next_rank, rank, MPI_COMM_WORLD, ierr) 
      call MPI_RECV(bufferx, size(bufferx), MPI_REAL, prev_rank, prev_rank, MPI_COMM_WORLD, status0, ierr)

      call MPI_SEND(sendy, size(bufferx), MPI_REAL, next_rank, rank, MPI_COMM_WORLD, ierr) 
      call MPI_RECV(buffery, size(bufferx), MPI_REAL, prev_rank, prev_rank, MPI_COMM_WORLD, status0, ierr)
      
      call MPI_SEND(sendz, size(bufferx), MPI_REAL, next_rank, rank, MPI_COMM_WORLD, ierr) 
      call MPI_RECV(bufferz, size(bufferx), MPI_REAL, prev_rank, prev_rank, MPI_COMM_WORLD, status0, ierr)
    end do

    deallocate(sendx)
    deallocate(sendy)
    deallocate(sendz)
  end function count_pairs
end program main
