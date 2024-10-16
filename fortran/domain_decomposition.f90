program main
  use iso_fortran_env
  use mpi

  implicit none

  ! mpi variables
  integer :: ierr, rank, nprocs

  ! cartesian communicator variables
  integer :: comm_cart
  integer :: ndims = 3
  integer, dimension(3) :: dims = [0, 0, 0]
  logical, dimension(3) :: periods = [.true.,.true.,.true.]
  logical :: reorder = .true.

  ! stores the position of the current process in the cartesian communicator grid
  integer, dimension(3) :: coord

  ! define boundaries and the cutoff distance
  real, dimension(3) :: lower_boundary = [0,0,0], upper_boundary = [1,1,1]
  real :: cutoff = 0.05

  ! arrays to store the particles positions and indexes
  real, dimension(:), allocatable :: posx, posy, posz
  integer, dimension(:), allocatable :: posi
  integer :: particle_count

  ! variables for domain decomposition
  real, dimension(3) :: lower_domain, upper_domain

  ! the number of pairs found
  integer(kind=int64) :: pairs, sum_pairs

  ! mpi boilerplate
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  ! setup timestamps
  if (rank == 0) print *, "[TIME] init", elapsed_time()

  ! find out how many domains per dimension
  call MPI_DIMS_CREATE(nprocs, ndims, dims, ierr)
  if (rank == 0) print *, "[TIME] dims created", elapsed_time() 
  if (rank == 0) print *, " [LOG] dims:", dims

  ! create the cartesian communicator and find the processes position in the grid
  call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm_cart, ierr)
  call MPI_CART_COORDS(comm_cart, rank, ndims, coord, ierr)
  if (rank == 0) print *, "[TIME] cartesian communicator setup", elapsed_time()

  ! find the domain bounds and filter the file for particles within the bounds
  call read_files(posx, posy, posz, posi, particle_count, &
                 lower_boundary, upper_boundary, lower_domain, upper_domain, cutoff, &
                 rank, dims, coord, comm_cart, ierr)
  if (rank == 0) print *, "[TIME] files read", elapsed_time()

  pairs = count_pairs( &
    posx, posy, posz, posi, particle_count, &
    lower_domain, upper_domain, &
    lower_boundary, upper_boundary, &
    cutoff, rank, coord, dims, comm_cart) 

  if (rank == 0) then
    print *, "[TIME] counted pairs", elapsed_time()
    print *, "[TIME] total elapsed time", elapsed_time(.true.)
  end if

  ! sum results
  call MPI_REDUCE(pairs, sum_pairs, 1, MPI_INTEGER8, MPI_SUM, 0, comm_cart, ierr)

  ! print results
  if (rank == 0) print *, " [LOG] pairs", sum_pairs

  ! end the program
  deallocate(posx)
  deallocate(posy)
  deallocate(posz)
  deallocate(posi)

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
  subroutine read_files(posx, posy, posz, posi, particle_count, &
    lower_boundary, upper_boundary, lower_domain, upper_domain, cutoff, &
    rank, dims, coord, comm, ierr)
    implicit none

    ! params
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer, dimension(:), allocatable, intent(inout) :: posi
    integer, intent(inout) :: particle_count

    real, dimension(:), allocatable :: tempx, tempy, tempz
    integer, dimension(:), allocatable :: tempi

    real, dimension(:), allocatable :: filterx, filtery, filterz
    integer, dimension(:), allocatable :: filteri

    real, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    real, dimension(3), intent(inout) :: lower_domain, upper_domain
    real, intent(inout) :: cutoff

    ! mpi variables
    integer, intent(in) :: rank, comm
    integer, intent(inout) :: ierr
    integer, dimension(3), intent(in) :: dims, coord
    
    ! temp variables
    integer :: num_particles, seed
    logical :: file_exists

    call read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed, rank, comm, ierr)
    call get_domain(dims, coord, lower_boundary, upper_boundary, lower_domain, upper_domain)

    inquire(file="particle_data.txt", exist=file_exists)

    if (file_exists) then
      call read_data(tempx, tempy, tempz, tempi, lower_boundary, upper_boundary, rank, comm, ierr)
    else
      call generate_data(tempx, tempy, tempz, tempi, seed, num_particles, lower_boundary, upper_boundary, rank, comm, ierr)
    end if

    ! find the particles within the domain and then save them into the local array
    call filter_particles(tempx, tempy, tempz, tempi, size(tempi), &
      lower_domain, upper_domain, lower_boundary, upper_boundary, &
      filterx, filtery, filterz, filteri)


    particle_count = size(filteri)
    allocate(posx(particle_count))
    allocate(posy(particle_count))
    allocate(posz(particle_count))
    allocate(posi(particle_count))

    posx(1:particle_count) = filterx
    posy(1:particle_count) = filtery
    posz(1:particle_count) = filterz
    posi(1:particle_count) = filteri
  end subroutine read_files

  ! reads the data from the config file
  subroutine read_config(lower_boundary, upper_boundary, cutoff, num_particles, seed, rank, comm, ierr)
    implicit none

    real, dimension(3), intent(inout) :: lower_boundary, upper_boundary
    real, intent(inout) :: cutoff
    integer, intent(inout) :: num_particles, seed 

    ! mpi variables
    integer, intent(in) :: rank, comm
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

    call MPI_BCAST(num_particles, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_BCAST(seed, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_BCAST(num_particles, 1, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(lower_boundary, 3, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(upper_boundary, 3, MPI_REAL, 0, comm, ierr)
  end subroutine read_config

  ! generates particle data
  subroutine generate_data(posx, posy, posz, posi, seed, num_particles, lower_boundary, upper_boundary, rank, comm, ierr)
    implicit none

    integer, intent(in) :: seed, num_particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer, dimension(:), allocatable, intent(inout) :: posi
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! array to generate the random particles into
    real, dimension(:,:), allocatable :: particles 

    ! distance between boundary walls on each axis
    real, dimension(3) :: boundary_diff
    
    ! looping varialbes
    integer :: axis, i

    ! mpi variables
    integer, intent(in) :: rank, comm
    integer, intent(inout) :: ierr

    ! random calls require an array of length 8
    integer, dimension(8) :: seed_array
    seed_array = seed

    ! allocate space in all arrays
    allocate(posx(num_particles))
    allocate(posy(num_particles))
    allocate(posz(num_particles))
    allocate(posi(num_particles))

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
        posi(i) = i
      end do
    end if

    call MPI_BCAST(posx, num_particles, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(posy, num_particles, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(posz, num_particles, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(posi, num_particles, MPI_INTEGER, 0, comm, ierr)
  end subroutine generate_data
  
  ! reads the file and extracts the particles within the processes domain
  subroutine read_data(posx, posy, posz, posi, &
                       lower_boundary, upper_boundary, &
                       rank, comm, ierr)
    implicit none

    ! variables to be filled with particles and the total number of particles
    real, dimension(:), allocatable, intent(inout) :: posx, posy, posz
    integer, dimension(:), allocatable, intent(inout) :: posi

    ! the bounds of the simulation
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! temporary variables that hold unfiltered data from the file
    integer :: file_length, line

    ! mpi variables
    integer, intent(in) :: rank, comm
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
    call MPI_BCAST(file_length, 1, MPI_INTEGER, 0, comm, ierr)

    ! allocate space into the arrays to store the particles
    allocate(posx(file_length))
    allocate(posy(file_length))
    allocate(posz(file_length))
    allocate(posi(file_length))

    if (rank == 0) then
      ! read the positions of the particles into the unfiltered array
      do line = 1, file_length
        read(11, 9) posx(line), posy(line), posz(line)
        posi(line) = line
      end do

      ! close file
      close(11)
    end if

    ! sync the values read from the file between the processes
    call MPI_BCAST(posx, file_length, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(posy, file_length, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(posz, file_length, MPI_REAL, 0, comm, ierr)
    call MPI_BCAST(posi, file_length, MPI_INTEGER, 0, comm, ierr)
  end subroutine read_data

  ! filters an array of particles for particles within a given region
  subroutine filter_particles(origx, origy, origz, origi, particle_count, & 
      lower_region, upper_region, lower_boundary, upper_boundary, &
      filterx, filtery, filterz, filteri)
    implicit none

    ! the particles to filter
    real, dimension(:), allocatable, intent(in) :: origx, origy, origz
    integer, dimension(:), allocatable, intent(in) :: origi
    integer, intent(in) :: particle_count

    ! the bounds of the region
    real, dimension(3), intent(in) :: lower_region, upper_region

    ! the bounds of the simulation
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    real, dimension(:), allocatable, intent(inout) :: filterx, filtery, filterz
    integer, dimension(:), allocatable, intent(inout) :: filteri

    ! temporary variables used for looping
    integer, dimension(:), allocatable :: indicies
    integer :: count, i
    
    ! allocate enough space for filtering all particles
    allocate(indicies(particle_count))

    ! find the particles that are withing the region and save their indicies
    count = 0
    do i = 1, particle_count
      if (within_region(origx(i), origy(i), origz(i), lower_region, upper_region)) then
        count = count + 1
        indicies(count) = i
      end if
    end do

    ! allocate the exact amount of space to store the filtered particles
    allocate(filterx(count))
    allocate(filtery(count))
    allocate(filterz(count))
    allocate(filteri(count))

    ! save the filtered particles into the arrays
    do i = 1, count
      filterx(i) = origx(indicies(i))
      filtery(i) = origy(indicies(i))
      filterz(i) = origz(indicies(i))
      filteri(i) = origi(indicies(i))
    end do
  end subroutine filter_particles

  ! calculates whether a position is within a domain. for each axis lower <= position < upper must be true
  logical function within_region(x, y, z, lower_region, upper_region) result(within)
    implicit none

    ! the coordinates the check
    real, intent(in) :: x, y, z

    ! the region the coordinate should be checked against
    real, dimension(3), intent(in) :: lower_region, upper_region

    ! variables for iteration
    real, dimension(3) :: position
    integer :: axis

    ! return true unless stated otherwise
    within = .true.
    
    ! store the coordinates in an indexed format
    position = [x, y, z]

    do axis = 1, 3
      if (position(axis) >= upper_region(axis) .or. position(axis) < lower_region(axis)) within = .false.
    end do
  end function within_region

  ! interpolates between two values, given a weight
  real function lerp(lower, upper, weight) result(res)
    implicit none
    real, intent(in) :: lower, upper, weight

    res = lower + weight * (upper - lower)
  end function lerp

  ! gets the domain of a process
  subroutine get_domain(dims, coord, lower_boundary, upper_boundary, lower_domain, upper_domain)
    implicit none

    integer, dimension(3), intent(in) :: dims, coord
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary
    real, dimension(3), intent(inout) :: lower_domain, upper_domain

    integer :: axis

    ! iterate over the axis and interpolate a distance between the upper and lower boundaries
    do axis = 1, size(dims)
      lower_domain(axis) = lerp(lower_boundary(axis), upper_boundary(axis), coord(axis) / real(dims(axis)))
      upper_domain(axis) = lerp(lower_boundary(axis), upper_boundary(axis), (coord(axis) + 1) / real(dims(axis)))
    end do
  end subroutine get_domain

  ! finds the distance between two points on an axis normally (ignoring PBC)
  real function distance(a, b) result(difference)
    implicit none

    ! a and b are two points along an axis
    real, intent(in) :: a, b

    ! return the difference between the two points
    difference = abs(a - b)
  end function distance

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

  ! checks whether two coordinates are a pair
  logical function check_pair(ax, ay, az, ai, bx, by, bz, bi, lower_boundary, upper_boundary, cutoff) result(is_pair)
    implicit none

    real, intent(in) :: ax, ay, az, bx, by, bz, cutoff
    integer, intent(in) :: ai, bi

    real, dimension(3), intent(in) :: lower_boundary, upper_boundary


    real, dimension(3) :: a, b, difference

    integer :: axis

    if (ai > bi) then
      ! store the values in a format that can be iterated over
      a = [ax, ay, az]
      b = [bx, by, bz]

      do axis = 1, 3
        difference(axis) = distance(a(axis), b(axis))
      end do

      ! check whether the magnitude of the difference is within the cutoff
      is_pair = magnitude(difference) < cutoff
    else
      is_pair = .false.
    end if
  end function check_pair

  subroutine correct_positions(direction, coord, dims, boundary_diff, particles)
    implicit none

    ! the direction the data is sent, the coord of the process and the number of dims on the axis
    integer, intent(in) :: direction, coord, dims

    ! the distance between the boundaries of the simulation on the axis
    real, intent(in) :: boundary_diff

    ! the array of particle positions on the axis
    real, dimension(:), intent(inout) :: particles
    
    ! looping variables
    integer :: i

    if (direction == -1 .and. coord == 0) then
      particles = particles - boundary_diff
    else if (direction == 1 .and. coord == dims - 1) then
      particles = particles + boundary_diff
    end if
  end subroutine correct_positions

  ! shares particles with neighbouring processes and counts the pairs between them
  integer(kind=int64) function count_pairs(posx, posy, posz, posi, particle_count, &
    lower_domain, upper_domain, lower_boundary, upper_boundary, cutoff, &
    rank, coord, dims, comm_cart) result(pairs)

    implicit none

    ! the positions of the particles within the processes domain
    real, dimension(:), allocatable, intent(in) :: posx, posy, posz
    integer, dimension(:), allocatable, intent(in) :: posi
    integer, intent(in) :: particle_count

    ! the bounds of the processes domain
    real, dimension(3), intent(in) :: lower_domain, upper_domain

    ! the boundaries of the simulation
    real, dimension(3), intent(in) :: lower_boundary, upper_boundary

    ! the maximum distance two particles can be from eachother before they are no longer considered a pair
    real, intent(in) :: cutoff

    ! mpi variables
    integer, intent(in) :: rank, comm_cart, coord(3)
    integer, dimension(3) :: dims
    integer :: status(MPI_STATUS_SIZE)
    
    ! looping variables
    integer :: x, y, z, axis, i, j
    integer, dimension(3) :: direction

    ! region to filter particles in
    real, dimension(3) :: lower_region, upper_region

    ! variables to store the filtered particles
    real, dimension(:), allocatable :: filterx, filtery, filterz, filter
    integer, dimension(:), allocatable :: filteri

    ! variables to store the accumulated neighbour particles
    real, dimension(:), allocatable :: accumx, accumy, accumz
    integer, dimension(:), allocatable :: accumi
    integer :: accum_count

    ! variables to temporarily store particles received from neighbours
    real, dimension(:), allocatable :: bufferx, buffery, bufferz
    integer, dimension(:), allocatable :: bufferi

    ! variables for sending data
    integer :: buffer_size, filter_size

    ! the ranks of the neighbours in the direction
    integer :: pos_neighbour, neg_neighbour

    ! the distance between the bounds of the simulation for each axis
    real, dimension(3) :: boundary_diff

    ! find the distance between the bounds on each axis
    do axis = 1, 3
      boundary_diff(axis) = abs(upper_boundary(axis) - lower_boundary(axis))
    end do

    ! allocate space to accumulate the particles from neighbours
    allocate(accumx(size(posi) * 3))
    allocate(accumy(size(posi) * 3))
    allocate(accumz(size(posi) * 3))
    allocate(accumi(size(posi) * 3))
    accum_count = 0

    ! share space around my processor
    j = 0
    do x = -1, 1
      do y = -1, 1
        do z = -1, 1
          j = j + 1
          if (x == 0 .and. y == 0 .and. z == 0) then
            do i = 1, particle_count
              accum_count = accum_count + 1
              accumx(accum_count) = posx(i) 
              accumy(accum_count) = posy(i) 
              accumz(accum_count) = posz(i) 
              accumi(accum_count) = posi(i) 
            end do
          else
            ! find the neighbours for the current iteration
            ! neg_neighbour will be the rank that receives my filtered particles
            ! pos_neighbour will be the rank that sends me its filtered particles
            call MPI_CART_RANK(comm_cart, [coord(1) + x, coord(2) + y, coord(3) + z], pos_neighbour, ierr)
            call MPI_CART_RANK(comm_cart, [coord(1) - x, coord(2) - y, coord(3) - z], neg_neighbour, ierr)

            ! calculate the region in which to look for particles
            direction = [x, y, z]
            do axis = 1, 3
              if (direction(axis) == 0) then
                lower_region(axis) = lower_domain(axis)
                upper_region(axis) = upper_domain(axis)
              else if (direction(axis) == 1) then
                lower_region(axis) = lower_domain(axis)
                upper_region(axis) = lower_domain(axis) + cutoff
              else
                lower_region(axis) = upper_domain(axis) - cutoff
                upper_region(axis) = upper_domain(axis)
              end if
            end do
            
            ! filter for particles within the region 
            call filter_particles( &
              posx, posy, posz, posi, particle_count, &
              lower_region, upper_region, &
              lower_boundary, upper_boundary, &
              filterx, filtery, filterz, filteri &
            )

            ! send the filtered particles to the neg_neighbour/ receive filtered particles from pos_neighbour
            filter_size = size(filteri)
            call MPI_SENDRECV(filter_size, 1, MPI_INTEGER, neg_neighbour, j + rank, &
                              buffer_size, 1, MPI_INTEGER, pos_neighbour, j + pos_neighbour, &
                              comm_cart, status, ierr)

            ! allocate memory into the receiving buffer arrays
            allocate(bufferx(buffer_size))
            allocate(buffery(buffer_size))
            allocate(bufferz(buffer_size))
            allocate(bufferi(buffer_size))

            ! send/receive particle data from neighbours
            call MPI_SENDRECV(filterx, filter_size, MPI_REAL, neg_neighbour, j + rank, &
                              bufferx, buffer_size, MPI_REAL, pos_neighbour, j + pos_neighbour, &
                              comm_cart, status, ierr)

            call MPI_SENDRECV(filtery, filter_size, MPI_REAL, neg_neighbour, j + rank, &
                              buffery, buffer_size, MPI_REAL, pos_neighbour, j + pos_neighbour, &
                              comm_cart, status, ierr)

            call MPI_SENDRECV(filterz, filter_size, MPI_REAL, neg_neighbour, j + rank, &
                              bufferz, buffer_size, MPI_REAL, pos_neighbour, j + pos_neighbour, &
                              comm_cart, status, ierr)

            call MPI_SENDRECV(filteri, filter_size, MPI_INTEGER, neg_neighbour, j + rank, &
                              bufferi, buffer_size, MPI_INTEGER, pos_neighbour, j + pos_neighbour, &
                              comm_cart, status, ierr)

            ! displace the received particles to be around my domain
            call correct_positions(x, coord(1), dims(1), boundary_diff(1), bufferx)
            call correct_positions(y, coord(2), dims(2), boundary_diff(2), buffery)
            call correct_positions(z, coord(3), dims(3), boundary_diff(3), bufferz)

            ! add the filtered particles to the accumulated list of neighbour particles
            do i = 1, buffer_size
              accum_count = accum_count + 1
              accumx(accum_count) = bufferx(i) 
              accumy(accum_count) = buffery(i) 
              accumz(accum_count) = bufferz(i) 
              accumi(accum_count) = bufferi(i) 
            end do

            ! reset the filter and buffer arrays
            deallocate(filterx, filtery, filterz, filteri)
            deallocate(bufferx, buffery, bufferz, bufferi)
          end if
        end do
      end do
    end do

    if (rank == 0) print *, "[TIME] shared particles", elapsed_time()

    ! loops between all possible pairs of particles to find which ones are in range
    pairs = 0
    do i = 1, particle_count 
      do j = 1, accum_count
        if (check_pair( &
          posx(i), posy(i), posz(i), posi(i), &
          accumx(j), accumy(j), accumz(j), accumi(j), &
          lower_boundary, upper_boundary, cutoff &
        )) then
          pairs = pairs + 1
        end if
      end do
    end do
  end function count_pairs
end program main
