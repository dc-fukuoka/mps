! http://li.mit.edu/Stuff/CNSE/Paper/Koshizuka96Oka.pdf

  ! NaN check
#if 0
#define myassert(val, word) \
  if (isnan(val)) then; \
     write(6,*) word," value:",val,__FILE__,":",__LINE__; \
#ifdef __INTEL_COMPILER
     call tracebackqq; \
#else
     call backtrace; \
     stop; \
  end if
#endif

module params
  implicit none

  ! mps parameters, variables
 
  integer,parameter::ndims = 3
  integer,parameter::dp = kind(1.0d0)

  real(dp)::dt = 3.0d-3 ! time slice

  integer,parameter::is_fluid=1,is_wall=2,is_dummy=3
  integer,dimension(ndims)::nparticles_fluid = (/8, 16,16/)
  integer,dimension(ndims)::nparticles_box   = (/32,16,16/)
  integer,dimension(ndims)::thick_wall  = (/2,2,2/)
  integer,dimension(ndims)::thick_dummy = (/2,2,2/)

  real(dp),parameter::rho = 1.0d3
  real(dp),parameter::nu  = 1.0d-6
  
  real(dp),parameter::dist_particles   = 7.5d-2
  real(dp),parameter::thre_nd          = 0.97d0
  real(dp),parameter::relax_pressure   = 0.2d0

  real(dp),parameter::coef_restitution = 2.0d-1 ! coefficient of restitution for rigid body collision
  real(dp),parameter::dist_col         = dist_particles*0.5d0 ! collision distance

  real(dp),parameter::re_nd    = dist_particles*2.1d0 ! effective radius for number density
  real(dp),parameter::re_grad  = dist_particles*2.1d0 ! effective radius for gradient
  real(dp),parameter::re_lap   = dist_particles*3.1d0 ! effective radius for laplacian

  real(dp),parameter::g = -9.8065d0 ! gravity constant

  integer::nparticles
  
  real(dp)::n0_nd
  real(dp)::n0_grad
  real(dp)::n0_lap
  real(dp)::lambda
  
  ! CG method parameters

  integer::iter_max=10000
  real(dp)::tol = 1.0d-10 ! CG method tolerance

  integer::tstep_max  = 200
  integer::freq_write = 200
  
  type particle
     real(dp),allocatable,dimension(:,:)::pos,vel,acc
     real(dp),allocatable,dimension(:)::ndensity,pressure,min_pressure,source
     integer,allocatable,dimension(:)::particle_type
     logical,allocatable,dimension(:)::is_surface
  end type particle

  type(particle)::particles
end module params

module mysubs
contains

  subroutine calc_nparticles
    use params
    implicit none
    integer::nparticles_wall_dummy
    nparticles_wall_dummy = (nparticles_box(1)+2*(thick_wall(1)+thick_dummy(1))) * &
                            (nparticles_box(2)+2*(thick_wall(2)+thick_dummy(2))) * &
                            (nparticles_box(3)+2*(thick_wall(3)+thick_dummy(3))) - &
                             nparticles_box(1)*nparticles_box(2)*(nparticles_box(3)+thick_wall(3)+thick_dummy(3)) ! upper surface is opened
    nparticles = nparticles_wall_dummy + nparticles_fluid(1)*nparticles_fluid(2)*nparticles_fluid(3)
    write(6,*) "total number of particles:",nparticles
    
  end subroutine calc_nparticles

  subroutine allocate_arrays
    use params
    implicit none
    allocate(particles%pos(nparticles,ndims))
    allocate(particles%vel(nparticles,ndims))
    allocate(particles%acc(nparticles,ndims))
    allocate(particles%ndensity(nparticles))
    allocate(particles%pressure(nparticles))
    allocate(particles%min_pressure(nparticles))
    allocate(particles%source(nparticles))
    allocate(particles%is_surface(nparticles))
    allocate(particles%particle_type(nparticles))

    !$omp parallel
    !$omp workshare
    particles%pos(:,:)         = 0.0d0
    particles%vel(:,:)         = 0.0d0
    particles%acc(:,:)         = 0.0d0
    particles%ndensity(:)      = 0.0d0
    particles%pressure(:)      = 0.0d0
    particles%min_pressure(:)  = 0.0d0
    particles%source(:)        = 0.0d0
    particles%is_surface(:)    = .false.
    particles%particle_type(:) = 0
    !$omp end workshare
    !$omp end parallel
    
  end subroutine allocate_arrays

  subroutine free_arrays
    use params
    implicit none
    deallocate(particles%pos)
    deallocate(particles%vel)
    deallocate(particles%acc)
    deallocate(particles%ndensity)
    deallocate(particles%pressure)
    deallocate(particles%min_pressure)
    deallocate(particles%source)
    deallocate(particles%is_surface)
    deallocate(particles%particle_type)
  end subroutine free_arrays

  subroutine set_initial_state
    use params
    implicit none
    integer::i,ix,iy,iz
    integer,dimension(ndims)::start,end
    integer::count
    real(dp)::x,y,z
    logical::is_particle=.false.
    integer::tot_dummy,tot_wall,tot_fluid,tot_empty
    
    do i=1,ndims
       start(i) = 1                 - thick_wall(i) - thick_dummy(i)
       end(i)   = nparticles_box(i) + thick_wall(i) + thick_dummy(i)
    end do

    count         = 1
    tot_dummy     = 0
    tot_wall      = 0
    tot_fluid     = 0
    tot_empty     = 0
    do iz=start(3),end(3)
       do iy=start(2),end(2)
          do ix=start(1),end(1)
             x = dist_particles*ix
             y = dist_particles*iy
             z = dist_particles*iz
             is_particle = .false.
             ! ok
             ! dummy region
             if (.not.( (ix.ge.start(1)+thick_dummy(1)).and.(ix.le.end(1)-thick_dummy(1))   .and. &
                        (iy.ge.start(2)+thick_dummy(2)).and.(iy.le.end(2)-thick_dummy(2))   .and. &
                        (iz.ge.start(3)+thick_dummy(3)).and.(iz.le.end(3)-thick_dummy(3)) ) .and. &
                 .not.( (ix.ge.start(1)+thick_dummy(1)).and.(ix.le.end(1)-thick_dummy(1))   .and. & ! upper surface is opened
                        (iy.ge.start(2)+thick_dummy(2)).and.(iy.le.end(2)-thick_dummy(2))   .and. &
                        (iz.ge.end(3)  -thick_dummy(3)).and.(iz.le.end(3)) ) ) then
                is_particle = .true.
                tot_dummy = tot_dummy + 1
                particles%particle_type(count) = is_dummy
             end if
             ! wall region
             if ( (ix.ge.start(1)+thick_dummy(1)).and.(ix.le.end(1)-thick_dummy(1)) .and. &
                  (iy.ge.start(2)+thick_dummy(2)).and.(iy.le.end(2)-thick_dummy(2)) .and. &
                  (iz.ge.start(3)+thick_dummy(3)).and.(iz.le.end(3)               ) .and. &
           .not.( (ix.ge.1).and.(ix.le.nparticles_box(1))   .and. &
                  (iy.ge.1).and.(iy.le.nparticles_box(2))   .and. &
                  (iz.ge.1).and.(iz.le.nparticles_box(3)) ) .and. &
           .not.( (ix.ge.1).and.(ix.le.nparticles_box(1))   .and. &
                  (iy.ge.1).and.(iy.le.nparticles_box(2))   .and. &
                  (iz.ge.nparticles_box(3)+1).and.(iz.le.end(3)) ) ) then
                is_particle = .true.
                tot_wall = tot_wall + 1
                particles%particle_type(count) = is_wall
             end if
             ! fluid particle region
             if ((ix.ge.1).and.(ix.le.nparticles_fluid(1)) .and. &
                 (iy.ge.1).and.(iy.le.nparticles_fluid(2)) .and. &
                 (iz.ge.1).and.(iz.le.nparticles_fluid(3))) then
                is_particle = .true.
                tot_fluid = tot_fluid + 1
                particles%particle_type(count) = is_fluid
             end if
             ! empty region
             if ( ( (ix.ge.1).and.(ix.le.nparticles_box(1))       .and. &
                    (iy.ge.1).and.(iy.le.nparticles_box(2))       .and. &
                    (iz.ge.1).and.(iz.le.nparticles_box(3))       .and. &
             .not.( (ix.ge.1).and.(ix.le.nparticles_fluid(1))     .and. &
                    (iy.ge.1).and.(iy.le.nparticles_fluid(2))     .and. &
                    (iz.ge.1).and.(iz.le.nparticles_fluid(3)) ) ) .or.  &
                  ( (ix.ge.1).and.(ix.le.nparticles_box(1))       .and. &
                    (iy.ge.1).and.(iy.le.nparticles_box(2))       .and. &
                    (iz.ge.nparticles_box(3)+1).and.(iz.le.end(3)) ) ) then
                is_particle = .false.
                tot_empty = tot_empty + 1
             end if

             if (is_particle) then
                particles%pos(count,1) = x
                particles%pos(count,2) = y
                particles%pos(count,3) = z
                particles%vel(count,1) = 0.0d0
                particles%vel(count,2) = 0.0d0
                particles%vel(count,3) = 0.0d0
                count = count + 1
             end if
          end do
       end do
    end do
    count = count - 1

    if (count.ne.tot_dummy+tot_wall+tot_fluid) then
       write(6,*) "count:",count
       write(6,*) "tot_dummy+tot_wall+tot_fluid:",tot_dummy+tot_wall+tot_fluid
#ifdef __INTEL_COMPILER
       call tracebackqq
#else
       call backtrace
#endif
       stop
    end if

    write(6,*) "# of dummy particles:",tot_dummy
    write(6,*) "# of wall particles: ",tot_wall
    write(6,*) "# of fluid particles:",tot_fluid
    
  end subroutine set_initial_state

  subroutine calc_n0_and_lambda
    use params
    implicit none
    integer::i
    integer::n_fluid1

    n0_nd   = 0.0d0
    n0_grad = 0.0d0
    n0_lap  = 0.0d0
    lambda  = 0.0d0

    n_fluid1 = find_fluid_start()
#ifdef _DEBUG    
    write(6,*) "n_fluid1:",n_fluid1
#endif

    ! need to choose a particle that is located almost of center of the fluid, not to have low density
    !omp parallel do reduction(+:n0_nd,n0_grad,n0_lap,lambda)
    do i=1,nparticles
       if (i.eq.n_fluid1) cycle
       n0_nd   = n0_nd   + weight(distance(i,n_fluid1),re_nd)
       n0_grad = n0_grad + weight(distance(i,n_fluid1),re_grad)
       n0_lap  = n0_lap  + weight(distance(i,n_fluid1),re_lap)
       lambda  = lambda  + (distance(i,n_fluid1)**2)*weight(distance(i,n_fluid1),re_lap)
    end do
    lambda = lambda/n0_lap
#ifdef _DEBUG
    write(6,"(a,4(1pe14.5))") "n0_nd,n0_grad,n0_lap,lambda:",n0_nd,n0_grad,n0_lap,lambda
#endif
  end subroutine calc_n0_and_lambda

  ! find the position of the first fluid particle
  ! it is located in dense position because of the corner
  function find_fluid_start() result(res)
    use params
    implicit none
    integer::res
    integer::i
    res   = 0
    do i=1,nparticles
       if (particles%particle_type(i).eq.is_fluid) then
          exit
       end if
    end do
    res = i
  end function find_fluid_start

  subroutine calc_ndensity
    use params
    implicit none
    integer::i,j

    !$omp parallel
    !$omp workshare
    particles%ndensity(:) = 0.0d0
    !$omp end workshare

    ! a better way is to find neighbor particles and only use them
    ! this calculation includes all of particles
    ! this includes all types of the particles
    !$omp do
    do i=1,nparticles
       do j=1,nparticles
          if (i.eq.j.or.distance(i,j).ge.re_nd) cycle
          particles%ndensity(i) = particles%ndensity(i) + weight(distance(i,j),re_nd)
       end do
    end do
    !$omp end do
    !$omp end parallel
    
  end subroutine calc_ndensity

  pure function weight(r,re) result(res)
    use params
    implicit none
    real(dp),intent(in)::r,re
    real(dp)::res

    if (r.ge.re) then
       res = 0.0d0
    else
       res = re/r - 1.0d0
    end if
    
  end function weight

  function distance(i,j) result(res)
    use params
    implicit none
    integer,intent(in)::i,j
    real(dp)::res

    res = sqrt((particles%pos(i,1)-particles%pos(j,1))**2 + &
               (particles%pos(i,2)-particles%pos(j,2))**2 + &
               (particles%pos(i,3)-particles%pos(j,3))**2)
  end function distance

  subroutine calc_gravity
    use params
    implicit none
    integer::i

    !$omp parallel
    !$omp workshare
    particles%acc(:,:) = 0.0d0
    !$omp end workshare
    !$omp do
    do i=1,nparticles
       if (particles%particle_type(i).eq.is_fluid) then
          particles%acc(i,3) = g ! z direction
       else
          particles%acc(i,:) = 0.0d0
       end if
    end do
    !$omp end do
    !$omp end parallel
  end subroutine calc_gravity

  subroutine calc_laplacian(in,out,type)
    use params
    implicit none
    real(dp),dimension(nparticles),intent(in)::in
    real(dp),dimension(nparticles),intent(out)::out
    character(1),intent(in)::type
    integer::i,j
    real(dp)::coef

    !$omp parallel
    !$omp workshare
    out(:) = 0.0d0
    !$omp end workshare
    !$omp end parallel

    coef = 2.0d0*ndims/n0_lap/lambda

    if (type.eq.'v') then ! for viscosity calculation
       !$omp parallel do
       do i=1,nparticles
          if (particles%particle_type(i).ne.is_fluid) cycle
          do j=1,nparticles
             if (i.eq.j.or.distance(j,i).gt.re_lap) cycle
             out(i) = out(i) + coef*(in(j)-in(i))*weight(distance(j,i),re_lap)
          end do
       end do
    else if (type.eq.'p') then ! for pressure calculation
       !$omp parallel do
       do i=1,nparticles
          ! presure is zero on the surface
          ! if (particles%is_surface(i)) cycle is mandatory for having positive definite matrix
          if (particles%is_surface(i).or.particles%particle_type(i).eq.is_dummy) cycle
          do j=1,nparticles
             if (i.eq.j.or.distance(j,i).ge.re_lap.or.particles%particle_type(j).eq.is_dummy) cycle
             out(i) = out(i) + coef*(in(j)-in(i))*weight(distance(j,i),re_lap)
          end do
       end do
    else
       write(6,*) "third argument of calc_laplacian() should be 'v' or 'p'."
       write(6,*) "third argument:",type
#ifdef __INTEL_COMPILER
       call tracebackqq
#else
       call backtrace
#endif
       stop
    end if
    
  end subroutine calc_laplacian

  ! calculate nu*laplacian(u)+g and substitute to acc(:,:)
  subroutine calc_viscosity
    use params
    implicit none
    integer::i,n
    real(dp),dimension(nparticles,ndims)::visc

    do n=1,ndims
       call calc_laplacian(particles%vel(1,n),visc(1,n),'v')
    end do
    
    !$omp parallel do
    do i=1,nparticles
       if (particles%particle_type(i).ne.is_fluid) cycle
       do n=1,ndims
          particles%acc(i,n) = particles%acc(i,n) + nu*visc(i,n)
       end do
    end do

  end subroutine calc_viscosity

  subroutine move_particles
    use params
    implicit none
    integer::i,n

    !$omp parallel do
    do i=1,nparticles
       do n=1,ndims
          if (particles%particle_type(i).ne.is_fluid) cycle
          particles%vel(i,n) = particles%vel(i,n) + particles%acc(i,n)*dt
          particles%pos(i,n) = particles%pos(i,n) + particles%vel(i,n)*dt
       end do
    end do

    !$omp parallel
    !$omp workshare
    particles%acc(:,:) = 0.0d0
    !$omp end workshare
    !$omp end parallel
    
  end subroutine move_particles

  subroutine calc_pressure
    use params
    implicit none

    call calc_ndensity
    call check_surface
    call calc_source
    call cg
    call remove_negative_pressure
    call calc_min_pressure
    
  end subroutine calc_pressure

  subroutine calc_min_pressure
    use params
    implicit none
    integer::i,j

    !$omp parallel do
    do i=1,nparticles
       if (particles%particle_type(i).eq.is_dummy) cycle
       particles%min_pressure(i) = particles%pressure(i)
       do j=1,nparticles
          if (i.eq.j.or.particles%particle_type(j).eq.is_dummy.or.distance(i,j).ge.re_grad) cycle
          if (particles%min_pressure(i).gt.particles%pressure(j)) then
             particles%min_pressure(i) = particles%pressure(j)
          end if
       end do
    end do
    
  end subroutine calc_min_pressure

  subroutine remove_negative_pressure
    use params
    implicit none
    integer::i

    !$omp parallel do
    do i=1,nparticles
       if (particles%pressure(i).lt.0.0d0) particles%pressure(i) = 0.0d0
    end do
  end subroutine remove_negative_pressure

  subroutine calc_gradient(in1,in2,out)
    use params
    implicit none
    real(dp),dimension(nparticles),intent(in)::in1,in2
    real(dp),dimension(nparticles,ndims),intent(out)::out
    real(dp)::coef
    integer::i,j,n

    !$omp parallel
    !$omp workshare
    out(:,:) = 0.0d0
    !$omp end workshare
    !$omp end parallel
    coef = ndims/n0_grad
    !$omp parallel do
    do i=1,nparticles
       if (particles%particle_type(i).ne.is_fluid) cycle
       do j=1,nparticles
          if (i.eq.j.or.particles%particle_type(i).eq.is_dummy.or.distance(j,i).ge.re_grad) cycle
          do n=1,ndims
             out(i,n) = out(i,n) + &
             coef*(in1(j)-in2(i))*(particles%pos(j,n)-particles%pos(i,n))/(distance(j,i)**2)*weight(distance(j,i),re_grad)
          end do
       end do
    end do
    
  end subroutine calc_gradient
  
  subroutine calc_a_dot_x(in,out)
    use params
    implicit none
    real(dp),dimension(nparticles),intent(in)::in
    real(dp),dimension(nparticles),intent(inout)::out
    integer::i

    !$omp parallel
    !$omp workshare
    out(:) = 0.0d0
    !$omp end workshare
    !$omp end parallel
    call calc_laplacian(in,out,'p')
    !$omp parallel do
    do i=1,nparticles
       out(i) = -1.0d0/rho*out(i)
    end do
    
  end subroutine calc_a_dot_x

  ! solve -1/rho*laplacian(p) = source for p
  subroutine cg
    use params
    implicit none
    real(dp),dimension(nparticles)::x,ax,x_new,p,ap,p_new,r,r_new,b
    real(dp)::alpha,beta,r2,r_new2,b2,eps,pap
    integer::i,iter

    !$omp parallel
    !$omp workshare
    b(:) = particles%source(:)
    x(:) = b(:) ! initial guess
    !$omp end workshare
    !$omp end parallel
    ! x(:) = 0.0d0

    call calc_a_dot_x(x,ax)
    !$omp parallel
    !$omp workshare
    r(:) = b(:) - ax(:)
    p(:) = r(:)
    !$omp end workshare
    !$omp end parallel

    b2 = 0.0d0
    !$omp parallel do reduction(+:b2)
    do i=1,nparticles
       b2 = b2 + b(i)**2
    end do

    do iter=1,iter_max
       alpha  = 0.0d0
       beta   = 0.0d0
       r2     = 0.0d0
       r_new2 = 0.0d0
       pap    = 0.0d0
       eps    = 0.0d0
       call calc_a_dot_x(p,ap)
       !$omp parallel do reduction(+:r2,pap)
       do i=1,nparticles
          r2  = r2  + r(i)**2
          pap = pap + p(i)*ap(i)
       end do
       if (pap.le.0.0d0) then
          write(6,*) "the matrix is not symmetric positive definite."
          write(6,*) "pap:",pap
#ifdef __INTEL_COMPILER
          call tracebackqq
#else
          call backtrace
#endif
          stop
       end if
       alpha = r2/pap
       !$omp parallel do reduction(+:r_new2)
       do i=1,nparticles
          x_new(i) = x(i)   + alpha*p(i)
          r_new(i) = r(i)   - alpha*ap(i)
          r_new2   = r_new2 + r_new(i)**2
       end do
       eps = sqrt(r_new2)/sqrt(b2)
       if (eps.le.tol) then
          exit
       else if (eps.gt.tol.and.iter.eq.iter_max) then
          write(6,*) "cg method did not converge."
          write(6,*) "eps,tol:",eps,tol
#ifdef __INTEL_COMPILER
          call tracebackqq
#else
          call backtrace
#endif
          stop
       else if (isnan(eps)) then
          write(6,*) "eps is NaN."
#ifdef __INTEL_COMPILER        
          call tracebackqq
#else
          call backtrace
#endif
          stop
       end if

       beta = r_new2/r2
       
       if (isnan(alpha).or.isnan(beta)) then
          write(6,*) "alpha or beta is NaN."
#ifdef __INTEL_COMPILER
          call tracebackqq
#else
          call backtrace
#endif
          stop
       end if

       !$omp parallel
       !$omp do
       do i=1,nparticles
          p_new(i) = r_new(i) + beta*p(i)
       end do
       !$omp end do
       !$omp do
       do i=1,nparticles
          x(i) = x_new(i)
          p(i) = p_new(i)
          r(i) = r_new(i)
       end do
       !$omp end do
       !$omp end parallel
#ifdef _DEBUG
       if (mod(iter,10).eq.0) then
          write(6,"(a,i6,1pe14.5)") "iter,eps:",iter,eps
       end if
#endif
    end do

    ! converged
    !$omp parallel
    !$omp workshare
    particles%pressure(:) = x_new(:)
    !$omp end workshare
    !$omp end parallel
    
  end subroutine cg

  subroutine check_surface
    use params
    implicit none
    integer::i

    !$omp parallel
    !$omp workshare
    particles%is_surface(:) = .false.
    !$omp end workshare
    !$omp do
    do i=1,nparticles
       if (particles%ndensity(i)/n0_nd.lt.thre_nd.and.particles%particle_type(i).ne.is_dummy) then
          particles%is_surface(i) = .true.
       end if
    end do
    !$omp end do
    !$omp end parallel

  end subroutine check_surface

  subroutine calc_source
    use params
    implicit none
    integer::i

    !$omp parallel
    !$omp workshare
    particles%source(:) = 0.0d0
    !$omp end workshare
    !$omp do
    do i=1,nparticles
       ! dummy is excluded here
       if (particles%particle_type(i).eq.is_dummy) cycle
       if (particles%is_surface(i)) then
          ! pressure in the surface is zero
          particles%source(i) = 0.0d0
       else
          particles%source(i) = relax_pressure*1.0d0/(dt**2)*(particles%ndensity(i)-n0_nd)/n0_nd
       end if
    end do
    !$omp end do
    !$omp end parallel
  end subroutine calc_source

  subroutine check_nan(tstep)
    use params
    implicit none
    integer,intent(in)::tstep
    integer::i,n

    do n=1,ndims
       do i=1,nparticles
          if (isnan(particles%pos(i,n))) then
             write(6,*) "tstep,i,n,pos:",tstep,i,n,particles%pos(i,n)
#ifdef __INTEL_COMPILER
             call tracebackqq
#else
             call backtrace
#endif
             stop
          else if (isnan(particles%vel(i,n))) then
             write(6,*) "tstep,i,n,vel:",tstep,i,n,particles%vel(i,n)
#ifdef __INTEL_COMPILER
             call tracebackqq
#else
             call backtrace
#endif
             stop
          end if
       end do
    end do
    
  end subroutine check_nan
  
  ! main routine
  subroutine mps
    use params
    implicit none
    integer::i,tstep
    real(dp)::t

    do tstep=1,tstep_max
#ifdef _DEBUG
       call check_nan(tstep)
#endif
       call calc_gravity
       call calc_viscosity
       call move_particles
       call collision
       call calc_pressure
       call move_particles_by_pressure
       if (mod(tstep,tstep_max/freq_write).eq.0) then
          t = tstep*dt
          write(6,*) "tstep,t:",tstep,t
          call write_output
       end if
    end do
    
  end subroutine mps

  ! rigid body collision if two particles are very close
  subroutine collision
    use params
    implicit none
    integer::i,j,n
    real(dp),dimension(nparticles,ndims)::vel_after_col
    real(dp),dimension(ndims)::vel_tmp
    real(dp)::fdt

    !$omp parallel
    !$omp do
    do i=1,nparticles
       do n=1,ndims
          vel_after_col(i,n) = particles%vel(i,n)
       end do
    end do
    !$omp end do
    !$omp do private(fdt,vel_tmp)
    do i=1,nparticles
       if (particles%particle_type(i).ne.is_fluid) cycle
       do n=1,ndims
          vel_tmp(n) = particles%vel(i,n)
       end do
       do j=1,nparticles
          if (i.eq.j.or.distance(i,j).ge.dist_col) cycle
          fdt = 0.0d0
          do n=1,ndims
             fdt = fdt + (particles%vel(i,n)-particles%vel(j,n))*(particles%pos(j,n)-particles%pos(i,n))/distance(i,j)
          end do
          if (fdt.gt.0.0d0) then
             fdt = fdt*(1.0d0+coef_restitution)*rho/2.0d0
             do n=1,ndims
                vel_tmp(n) = vel_tmp(n) - fdt/rho*(particles%pos(j,n)-particles%pos(i,n))/distance(i,j)
             end do
          end if
#ifdef _DEBUG
          if (j.gt.i) then
             write(6,"(a,i0,a,i0,a,1pe14.5)") "collision occurred between ",i," and ",j,", distance:",distance(i,j)
          end if
#endif
       end do
       do n=1,ndims
          vel_after_col(i,n) = vel_tmp(n)
       end do
    end do
    !$omp end do
    !$omp do
    do i=1,nparticles
       if (particles%particle_type(i).ne.is_fluid) cycle
       do n=1,ndims
          particles%pos(i,n) = particles%pos(i,n)+(vel_after_col(i,n)-particles%vel(i,n))*dt
          particles%vel(i,n) = vel_after_col(i,n)
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine collision

  subroutine write_output
    use params
    implicit none
    integer::i
    
    do i=1,nparticles
       if (particles%particle_type(i).ne.is_fluid) cycle
       write(999,"(3(1pe14.5))") particles%pos(i,1),particles%pos(i,2),particles%pos(i,3)
    end do
    write(999,*)
    write(999,*)

  end subroutine write_output

  subroutine move_particles_by_pressure
   use params
   implicit none
   real(dp),dimension(nparticles,ndims)::grad_pres
   integer::i,n

   !$omp parallel
   !$omp workshare
   grad_pres(:,:) = 0.0d0
   !$omp end workshare
   !$omp end parallel
   ! corrected gradient model, use minimum pressure in the effective radius
   ! this guarantees threre is no negative pressure
   call calc_gradient(particles%pressure,particles%min_pressure,grad_pres)

   !$omp parallel
   !$omp do
   do i=1,nparticles
      if (particles%particle_type(i).ne.is_fluid) cycle
      do n=1,ndims
         particles%acc(i,n) = -1.0d0/rho*grad_pres(i,n)
      end do
   end do
   !$omp end do
#if 0
#ifdef _DEBUG
   do i=1,nparticles
      if (particles%particle_type(i).ne.is_fluid) cycle
      write(888,"(l2,3(1pe14.5))") particles%is_surface(i),(particles%acc(i,n),n=1,3)
   end do
#endif
#endif
   !$omp do
   do i=1,nparticles
      do n=1,ndims
         if (particles%particle_type(i).ne.is_fluid) cycle
         particles%vel(i,n) = particles%vel(i,n) + particles%acc(i,n)*dt
         particles%pos(i,n) = particles%pos(i,n) + particles%acc(i,n)*(dt**2)
      end do
   end do
   !$omp end do
   !$omp workshare
   particles%acc(:,:) = 0.0d0
   !$omp end workshare
   !$omp end parallel
   
 end subroutine move_particles_by_pressure
  
  subroutine debug
    use params
    implicit none
    integer::i,count_fluid,count_else

    do i=1,nparticles
       if (particles%particle_type(i).eq.is_fluid) then
          write(100,"(3(1pe14.5))") particles%pos(i,1), particles%pos(i,2), particles%pos(i,3)
          count_fluid = count_fluid + 1
       else if (particles%particle_type(i).eq.is_dummy.or.particles%particle_type(i).eq.is_wall) then
          write(101,"(3(1pe14.5))") particles%pos(i,1), particles%pos(i,2), particles%pos(i,3)
          count_else = count_else + 1
       end if
    end do
    
    write(6,*) "count_fluid,count_else:",count_fluid,count_else
    
  end subroutine debug

#ifndef __INTEL_COMPILER
  function dclock() result(res)
    use params
    implicit none
    real(dp)::res,dble
    integer::clk,clk_rate,clk_max
    
    call system_clock(clk,clk_rate,clk_max)
    res = dble(clk)/dble(clk_rate)
  end function dclock
#endif
  
end module mysubs

program main
  use params
  use mysubs
  implicit none
#ifdef __INTEL_COMPILER
  real(dp)::dclock
#endif
  real(dp)::time,t0

  call calc_nparticles
  call allocate_arrays
  call set_initial_state
  call calc_n0_and_lambda
  t0 = dclock()
  call mps
  time = dclock() - t0
  write(6,*) "time[s]:", time
!  call debug
  call free_arrays
  
  stop
end program main
