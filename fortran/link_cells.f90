program main
  use iso_fortran_env
  implicit none

  ! define positions of the particles
  double precision, dimension(:), allocatable :: posx, posy, posz

  ! define arrays for constructing link cells
  integer, dimension(:,:,:), allocatable :: head_of_chain
  integer, dimension(:), allocatable :: next_index
  double precision, dimension(3) :: cell_size
  integer, dimension(3) :: boundary_divisions

  ! define boundaries of the simulation and the cutoff distance
  double precision, dimension(3) :: lower_boundary, upper_boundary
  double precision :: cutoff

  ! stores the number of pairs counted
  integer(kind=int64) :: pairs

  print *, "[TIME] init", elapsed_time()

  ! read the file
  call read_files(posx, posy, posz, lower_boundary, upper_boundary, cutoff)
  print *, "[TIME] files read", elapsed_time()

  call get_cell_size(lower_boundary, upper_boundary, cutoff, cell_size, boundary_divisions) 
  print *, " [LOG] cell size", cell_size
  print *, " [LOG] boundary divisions", boundary_divisions

  call construct_link_cells( &
    posx, posy, posz, &
    cell_size, boundary_divisions, &
    lower_boundary, upper_boundary, &
    head_of_chain, next_index)

  ! run the pair finding algorithm and print the result
  pairs = count_pairs(posx, posy, posz, cutoff, head_of_chain, next_index, boundary_divisions)
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

  subroutine get_cell_size(lower_boundary, upper_boundary, cutoff, cell_size, boundary_divisions)
    implicit none

    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary
    double precision, intent(in) :: cutoff
    double precision, dimension(3), intent(inout) :: cell_size
    integer, dimension(3), intent(inout) :: boundary_divisions

    integer :: axis

    do axis = 1, 3
      boundary_divisions(axis) = int(abs(upper_boundary(axis) - lower_boundary(axis)) / cutoff)
      cell_size(axis) = abs(upper_boundary(axis) - lower_boundary(axis)) / dble(boundary_divisions(axis))
    end do
  end subroutine get_cell_size

  subroutine construct_link_cells( &
    posx, posy, posz, &
    cell_size, boundary_divisions, &
    lower_boundary, upper_boundary, &
    head_of_chain, next_index)
    implicit none

    double precision, dimension(:), allocatable, intent(in) :: posx, posy, posz

    double precision, dimension(3), intent(in) :: cell_size
    integer, dimension(3), intent(in) :: boundary_divisions

    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary

    integer, dimension(:,:,:), allocatable, intent(inout) :: head_of_chain
    integer, dimension(:), allocatable, intent(inout) :: next_index

    integer, dimension(3) :: cell_coord
    integer :: cell_index, i, j

    allocate(head_of_chain(0:boundary_divisions(1)-1, 0:boundary_divisions(2)-1, 0:boundary_divisions(3)-1))
    allocate(next_index(size(posx)))

    head_of_chain = 0
    next_index = 0

    do i = 1, size(posx)
      call get_cell_coord([posx(i), posy(i), posz(i)], boundary_divisions, lower_boundary, upper_boundary, cell_coord)

      j = head_of_chain(cell_coord(1), cell_coord(2), cell_coord(3))
      if (j == 0) then
        head_of_chain(cell_coord(1), cell_coord(2), cell_coord(3)) = i
      else
        do
          if (next_index(j) == 0) then
            exit
          else
            j = next_index(j)
          end if
        end do
        next_index(j) = i
      end if
    end do
  end subroutine construct_link_cells

  subroutine get_cell_coord(position, boundary_divisions, lower_boundary, upper_boundary, coord)
    implicit none

    double precision, dimension(3), intent(in) :: position
    integer, dimension(3), intent(in) :: boundary_divisions
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary

    integer, dimension(3), intent(inout) :: coord

    integer :: axis

    do axis = 1, 3
      coord(axis) = floor(boundary_divisions(axis) * abs(position(axis) - lower_boundary(axis)) / abs(upper_boundary(axis) - lower_boundary(axis)))
    end do
  end subroutine get_cell_coord

  ! reads data from the config and particle data files
  subroutine read_files(posx, posy, posz, lower_boundary, upper_boundary, cutoff)
    implicit none

    double precision, dimension(:), allocatable, intent(inout) :: posx, posy, posz

    double precision, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    double precision, intent(inout) :: cutoff
    
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

    double precision, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    double precision, intent(inout) :: cutoff
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
    double precision, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    double precision, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! array to generate the random particles into
    double precision, dimension(:,:), allocatable :: particles 

    ! distance between boundary walls on each axis
    double precision, dimension(3) :: boundary_diff
    
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
    double precision, dimension(:), allocatable, intent(inout) :: posx, posy, posz
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
  double precision function squared_magnitude(vector) result(res)
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
  end function squared_magnitude

  integer(kind=int64) function count_pairs(posx, posy, posz, cutoff, head_of_chain, next_index, boundary_divisions) result(pairs)
    implicit none

    ! the positions of the particles and the boundaries they reside in
    double precision, dimension(:), allocatable, intent(in) :: posx, posy, posz
    integer, dimension(:,:,:), allocatable, intent(in) :: head_of_chain
    integer, dimension(:), allocatable, intent(in) :: next_index

    ! the maximum distance between particles that can be considered a "pair"
    double precision, intent(in) :: cutoff
    
    integer, dimension(3), intent(in) :: boundary_divisions

    ! temporary variables for looping
    integer :: di, i, j, axis, k

    integer :: current_x, current_y, current_z, current_index
    integer :: neighbour_x, neighbour_y, neighbour_z, neighbour_index

    ! current is the current particle in the loop
    double precision, dimension(3) :: i_pos, j_pos, difference

    ! the offsets to look at for link cells
    integer, dimension(13) :: dx, dy, dz, visited_neighbours

    dx = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 ]
    dy = [-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0 ]
    dz = [-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1, 1 ]

    pairs = 0
    do current_x = 0, boundary_divisions(1) - 1 
      do current_y = 0, boundary_divisions(2) - 1 
        do current_z = 0, boundary_divisions(3) - 1 
          i = head_of_chain(current_x, current_y, current_z)
          if (i == 0) cycle

          current_index = current_x + &
            current_y * boundary_divisions(1) + &
            current_z * boundary_divisions(1) * boundary_divisions(2)

          do
            j = next_index(i)
            if (j == 0) exit

            i_pos = [posx(i), posy(i), posz(i)]

            do
              j_pos = [posx(j), posy(j), posz(j)]
              
              ! find the difference between the two positions
              do axis = 1, 3
                difference(axis) = shortest_distance(i_pos(axis), j_pos(axis), lower_boundary(axis), upper_boundary(axis))
              end do

              ! if the difference between positions is within the cutoff count the pair
              if (squared_magnitude(difference) < cutoff ** 2) then
                pairs = pairs + 1
              end if

              if (next_index(j) == 0) then
                exit
              else
                j = next_index(j)
              end if
            end do

            i = next_index(i)
          end do

          visited_neighbours = 0
          neighbour_loop: do di = 1, size(dx)
            neighbour_x = mod(current_x + dx(di), boundary_divisions(1))
            if (neighbour_x < 0) neighbour_x = boundary_divisions(1) + neighbour_x

            neighbour_y = mod(current_y + dy(di), boundary_divisions(2))
            if (neighbour_y < 0) neighbour_y = boundary_divisions(2) + neighbour_y

            neighbour_z = mod(current_z + dz(di), boundary_divisions(3))
            if (neighbour_z < 0) neighbour_z = boundary_divisions(3) + neighbour_z

            neighbour_index = neighbour_x + &
              neighbour_y * boundary_divisions(1) + &
              neighbour_z * boundary_divisions(1) * boundary_divisions(2)
            visited_neighbours(di) = neighbour_index

            do k = 1, di - 1
              if (neighbour_index == visited_neighbours(k)) cycle neighbour_loop
            end do

            if (head_of_chain(neighbour_x, neighbour_y, neighbour_z) == 0) cycle

            i = head_of_chain(current_x, current_y, current_z)
            do
              i_pos = [posx(i), posy(i), posz(i)]

              j = head_of_chain(neighbour_x, neighbour_y, neighbour_z)
              do
                j_pos = [posx(j), posy(j), posz(j)]

                ! find the difference between the two positions
                do axis = 1, 3
                  difference(axis) = shortest_distance(i_pos(axis), j_pos(axis), lower_boundary(axis), upper_boundary(axis))
                end do

                ! if the difference between positions is within the cutoff count the pair
                if (squared_magnitude(difference) < cutoff ** 2) then
                  pairs = pairs + 1
                end if

                if (next_index(j) == 0) then
                  exit
                else
                  j = next_index(j)
                end if
              end do    

              if (next_index(i) == 0) then
                exit
              else
                i = next_index(i)
              end if
            end do
          end do neighbour_loop
        end do
      end do
    end do
  end function count_pairs
end program main
