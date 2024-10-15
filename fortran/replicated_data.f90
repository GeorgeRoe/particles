program main
  use iso_fortran_env
  use mpi

  implicit none

  ! mpi variables
  integer :: ierr, rank, nprocs

  ! define positions of the particles
  real, dimension(:), allocatable :: posx, posy, posz
  integer :: particle_count

  ! define boundaries of the simulation and the cutoff distance
  real, dimension(3) :: lower_boundary = [0,0,0], upper_boundary = [1,1,1]
  real :: cutoff = 0.05

  ! variables for storing the indicies in which the process should chech
  integer(kind=int64) :: start_index, end_index, total_pairs

  ! stores the number of pairs counted
  integer(kind=int64) :: pairs, sum_pairs

  ! mpi boilerplate
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  ! read the necessary files
  call read_files(posx, posy, posz, total_particle_count, lower_boundary, upper_boundary, cutoff, rank, ierr)

  ! find the starting and ending index that the process should work between
  ! calculations are split to avoid implicitly lowering precision
  total_pairs = particle_count - 1
  total_pairs = total_pairs * particle_count
  total_pairs = total_pairs / 2

  if (rank == 0) then
    start_index = 0
  else
    start_index = (rank / dble(nprocs)) * total_pairs + 1
  end if

  if (rank == nprocs - 1) then
    end_index = total_pairs
  else
    end_index = ((rank + 1) / dble(nprocs)) * total_pairs
  end if
  print *, "rank", rank, "will start at index", start_index, "and end at index", end_index, "of", total_pairs

  ! run the pair finding algorithm and print the result
  pairs = count_pairs( &
    posx, posy, posz, &
    lower_boundary, upper_boundary, &
    start_index, end_index, cutoff)
  print *, "rank", rank, "found", pairs, "pairs"

  call MPI_REDUCE(pairs, sum_pairs, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (rank == 0) print *, "total pairs found", sum_pairs

  deallocate(posx)
  deallocate(posy)
  deallocate(posz)

  call MPI_FINALIZE(ierr)

contains

  ! reads data from the config and particle data files
  subroutine read_files(posx, posy, posz, &
    particle_count, lower_boundary, upper_boundary, cutoff, &
    rank, ierr)
    implicit none

    ! params
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer, intent(inout) :: particle_count

    real, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    real, intent(inout) :: cutoff

    ! mpi variables
    integer, intent(in) :: rank
    integer, intent(inout) :: ierr
    
    ! temp variables
    integer :: num_particles, seed
    logical :: file_exists

    call read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed, rank, ierr)

    inquire(file="particle_data.txt", exist=file_exists)

    if (file_exists) then
      call read_data(posx, posy, posz, lower_boundary, upper_boundary, rank, ierr)
    else
      call generate_data(posx, posy, posz, seed, num_particles, lower_boundary, upper_boundary, rank, ierr)
    end if
  end subroutine read_files

  ! reads the data from the config file
  subroutine read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed, rank, ierr)
    implicit none

    real, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    real, intent(inout) :: cutoff
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
    call MPI_BCAST(num_particles, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(lower_boundary, 3, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(upper_boundary, 3, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_config

  ! generates particle data
  subroutine generate_data(posx, posy, posz, seed, num_particles, lower_boundary, upper_boundary, rank, ierr)
    implicit none

    integer, intent(in) :: seed, num_particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! array to generate the random particles into
    real, dimension(:,:), allocatable :: particles 

    ! distance between boundary walls on each axis
    real, dimension(3) :: boundary_diff
    
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
    end if

    call MPI_BCAST(posx, num_particles, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posy, num_particles, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posz, num_particles, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  end subroutine generate_data
  
  ! reads the file and extracts the particles within the processes domain
  subroutine read_data(posx, posy, posz, &
                       lower_boundary, upper_boundary, &
                       rank, ierr)
    implicit none

    ! variables to be filled with particles and the total number of particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz

    ! the bounds of the simulation
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

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
    call MPI_BCAST(posx, file_length, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posy, file_length, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(posz, file_length, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  end subroutine read_data

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

  integer(kind=int64) function count_pairs( &
      posx, posy, posz, &
      lower_boundary, upper_boundary, &
      start_index, end_index, cutoff &
      ) result(pairs)
    implicit none

    ! the positions of the particles and the boundaries they reside in
    real, dimension(:), allocatable, intent(in) :: posx, posy, posz
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! the pairs to iterate over and the indexes to start and end at
    integer(kind=int64), intent(in) :: start_index, end_index

    ! the maximum distance between particles that can be considered a "pair"
    real, intent(in) :: cutoff

    ! temporary variables for looping
    integer :: i, j, axis
    integer(kind=int64) :: current_pair

    ! current is the current particle in the loop
    real, dimension(3) :: i_pos, j_pos, difference

    pairs = 0
    current_pair = 0

    ! loop through the pairs
    do i = 1, size(posx)
      i_pos = [posx(i), posy(i), posz(i)]
      do j = i + 1, size(posx)
        ! if the current pair index is out of the bounds for the process ignore the pair
        current_pair = current_pair + 1
        if (current_pair >= start_index .and. current_pair <= end_index) then
          j_pos = [posx(j), posy(j), posz(j)]

          ! find the difference between the two positions
          do axis = 1, 3
            difference(axis) = shortest_distance( &
              i_pos(axis), j_pos(axis), &
              lower_boundary(axis), upper_boundary(axis))
          end do
          
          ! if the difference between positions is within the cutoff count the pair
          if (magnitude(difference) < cutoff) then
            pairs = pairs + 1
          end if
        end if
      end do
    end do
  end function count_pairs
end program main
