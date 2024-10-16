program main
  use iso_fortran_env
  implicit none

  ! define positions of the particles
  real, dimension(:), allocatable :: posx, posy, posz

  ! define boundaries of the simulation and the cutoff distance
  real, dimension(3) :: lower_boundary, upper_boundary
  real :: cutoff

  ! stores the number of pairs counted
  integer(kind=int64) :: pairs

  print *, "[TIME] init", elapsed_time()

  ! read the file
  call read_files(posx, posy, posz, lower_boundary, upper_boundary, cutoff)
  print *, "[TIME] files read", elapsed_time()

  ! run the pair finding algorithm and print the result
  pairs = count_pairs(posx, posy, posz, lower_boundary, upper_boundary, cutoff)
  print *, "[TIME] counted pairs", elapsed_time()  
  
  print *, "[TIME] total elapsed time", elapsed_time(.true.)
  print *, " [LOG] pairs", pairs

contains

  ! returns the elapsed time since the function was last called
  real function elapsed_time(get_total) result(elapsed)
    implicit none

    ! whether or not to return the running total
    logical, intent(in), optional :: get_total
    
    ! implicit static keeps values between function calls
    real :: start_time = 0.0, end_time = 0.0, running_total = 0.0
    logical :: initialised = .false.
  
    if (present(get_total) .and. initialised) then ! if asking for the total then return it
      elapsed = running_total
    else if (initialised) then ! if the function has already been called once before
      start_time = end_time
      call cpu_time(end_time)
      elapsed = end_time - start_time
      running_total = running_total + elapsed
    else ! if this is the first time the function is ran
      call cpu_time(start_time)
      initialised = .true.
      elapsed = 0.0
    end if
  end function elapsed_time

  ! reads data from the config and particle data files
  subroutine read_files(posx, posy, posz, lower_boundary, upper_boundary, cutoff)
    implicit none

    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz

    real, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    real, intent(inout) :: cutoff
    
    integer :: num_particles, seed

    logical :: file_exists

    call read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed)

    inquire(file="particle_data.txt", exist=file_exists)
    if (file_exists) then
      call read_data(posx, posy, posz)
    else
      call generate_data(seed, num_particles, lower_boundary, upper_boundary, posx, posy, posz)
    end if

  end subroutine read_files

  ! reads the data from the config file
  subroutine read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed)
    implicit none

    real, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    real, intent(inout) :: cutoff
    integer, intent(inout) :: num_particles, seed 

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
  end subroutine read_config

  ! generates particle data
  subroutine generate_data(seed, num_particles, lower_boundary, upper_boundary, posx, posy, posz)
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

    ! random calls require an array of length 8
    integer, dimension(8) :: seed_array
    seed_array = seed

    ! find the distance between the simulation bounds
    do axis = 1, 3
      boundary_diff(axis) = upper_boundary(axis) - lower_boundary(axis)
    end do

    ! allocate space in all arrays
    allocate(posx(num_particles))
    allocate(posy(num_particles))
    allocate(posz(num_particles))

    allocate(particles(num_particles,3))

    ! set the seed and then generate the random particles
    call random_seed(put=seed_array)
    call random_number(particles)

    do i = 1, num_particles
      posx(i) = lower_boundary(1) + particles(i, 1) * boundary_diff(1)
      posy(i) = lower_boundary(2) + particles(i, 2) * boundary_diff(2)
      posz(i) = lower_boundary(3) + particles(i, 3) * boundary_diff(3)
    end do
  end subroutine generate_data

  ! reads particle data from the particle data file
  subroutine read_data(posx, posy, posz)
    implicit none

    ! files are formatted as such:
    ! the first line contains how many rows worth of data the file stores
    ! the subsequent lines store data in columns for the x, y and z positions of the particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer :: file_length, line

    ! open the file
    open(11, file = "particle_data.txt", status = "old")

    ! read how many lines of data the file holds and allocate the space
    read(11, *) file_length

    allocate(posx(file_length))
    allocate(posy(file_length))
    allocate(posz(file_length))

    ! define the format of the columnar data and read it into the variables
    9 format(f20.17, 4x, f20.17, 4x, f20.17)
    do line = 1, file_length
      read(11, 9) posx(line), posy(line), posz(line)
    end do

    close(11)
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

  integer(kind=int64) function count_pairs(posx, posy, posz, lower_boundary, upper_boundary, cutoff) result(pairs)
    implicit none

    ! the positions of the particles and the boundaries they reside in
    real, dimension(:), allocatable, intent(in) :: posx, posy, posz
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! the maximum distance between particles that can be considered a "pair"
    real, intent(in) :: cutoff

    ! temporary variables for looping
    integer :: i, j, axis

    ! current is the current particle in the loop
    real, dimension(3) :: i_pos, j_pos, difference

    pairs = 0
  
    ! nested loop to get each combination of particles
    do i = 1, size(posx)
      ! get the current positions of the particles being looped over
      i_pos = [posx(i), posy(i), posz(i)]
      do j = i + 1, size(posx)
        j_pos = [posx(j), posy(j), posz(j)]

        ! find the difference between the two positions
        do axis = 1, 3
          difference(axis) = shortest_distance(i_pos(axis), j_pos(axis), lower_boundary(axis), upper_boundary(axis))
        end do

        ! if the difference between positions is within the cutoff count the pair
        if (magnitude(difference) < cutoff) then
          pairs = pairs + 1
        end if
      end do
    end do
  end function count_pairs
end program main
