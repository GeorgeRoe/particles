program main
  use iso_fortran_env
  implicit none

  ! define positions of the particles
  real, dimension(:), allocatable :: posx, posy, posz

  ! define boundaries of the simulation and the cutoff distance
  real, dimension(3) :: lower_boundary = [0,0,0], upper_boundary = [1,1,1]
  real :: cutoff = 0.5

  ! stores the number of pairs counted
  integer(kind=int64) :: pairs

  ! read the file
  call read_file(posx, posy, posz)

  ! run the pair finding algorithm and print the result
  pairs = count_pairs(posx, posy, posz, lower_boundary, upper_boundary, cutoff)

  print *, pairs

contains
  ! reads particle data from the particle data file
  subroutine read_file(posx, posy, posz)
    implicit none

    ! files are formatted as such:
    ! the first line contains how many rows worth of data the file stores
    ! the subsequent lines store data in columns for the x, y and z positions of the particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer :: file_length

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
