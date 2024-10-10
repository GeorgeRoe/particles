program main
	implicit none

	integer :: numparticles, counter, seed_input

	integer, dimension(8) :: seed
	integer :: nParticles
	real, dimension(3) :: axisDimension = [1.0, 1.0, 1.0]
	real, dimension(3) :: minBoundary = [0.0, 0.0, 0.0]

	real, dimension(:,:), allocatable :: Aparticles
	real, dimension(:,:), allocatable :: particles

	print *, "Number of particles:"
	read(*,*) nParticles

	print *, "Seed:"
	read(*,*) seed_input

	seed = seed_input

	allocate(Aparticles(nParticles,3))
	allocate(particles(nParticles,3))

	call random_seed(put=seed)
	call random_number(Aparticles)
	
	do counter = 1, nParticles
		particles(counter, 1:3) = [&
			& minBoundary(1) + Aparticles(counter, 1)*axisDimension(1), &
			& minBoundary(2) + Aparticles(counter, 2)*axisDimension(2), &
			& minBoundary(3) + Aparticles(counter, 3)*axisDimension(3) &
		&]
	end do


	open(10, file="particle_data.txt", status="old")

	9 format(f20.17, 4x, f20.17, 4x, f20.17)

	write(10, *) nParticles

	do counter = 1, nParticles
		write(10, 9) particles(counter,1), particles(counter,2), particles(counter,3)
	end do

	close(10)

	print *, "Generated"
	deallocate(Aparticles, particles)
end program main
