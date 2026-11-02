! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0

module test_hessian
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use multicharge, only: get_charges, add_charge_hessian, cache_container
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_eeq, only: eeq_cache
   use multicharge_model_hessian, only: coulomb_kernel, radial_hessian
   use multicharge_model_type, only: get_dir_trans
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private
   public :: collect_hessian
contains

subroutine collect_hessian(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      & new_unittest("eeq-molecule", test_eeq), &
      & new_unittest("eeqbc-molecule", test_eeqbc), &
      & new_unittest("eeq-cn-cap", test_eeq_cap), &
      & new_unittest("eeqbc-cn-cap", test_eeqbc_cap), &
      & new_unittest("eeq-charged", test_eeq_charged), &
      & new_unittest("eeqbc-charged", test_eeqbc_charged), &
      & new_unittest("eeq-periodic", test_eeq_periodic), &
      & new_unittest("eeqbc-periodic", test_eeqbc_periodic), &
      & new_unittest("isolated-atom", test_atom), &
      & new_unittest("response-cache", test_cache), &
      & new_unittest("invalid-response", test_invalid), &
      & new_unittest("small-r-kernel-limits", test_small_r_kernel) &
      & ]
end subroutine collect_hessian

subroutine test_eeq(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 1)
end subroutine test_eeq

subroutine test_eeqbc(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .true., 1)
end subroutine test_eeqbc

subroutine test_eeq_cap(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 2)
end subroutine test_eeq_cap

subroutine test_eeqbc_cap(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .true., 2)
end subroutine test_eeqbc_cap

subroutine test_eeq_charged(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 3)
end subroutine test_eeq_charged

subroutine test_eeqbc_charged(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .true., 3)
end subroutine test_eeqbc_charged

subroutine test_eeq_periodic(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: icase
   do icase = 5, 9
      call run_case(error, .false., icase)
      if (allocated(error)) return
   end do
   do icase = 12, 15
      call run_case(error, .false., icase)
      if (allocated(error)) return
   end do
end subroutine test_eeq_periodic

subroutine test_eeqbc_periodic(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: icase
   do icase = 5, 15
      call run_case(error, .true., icase)
      if (allocated(error)) return
   end do
end subroutine test_eeqbc_periodic

subroutine test_atom(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 4)
   if (allocated(error)) return
   call run_case(error, .true., 4)
end subroutine test_atom

subroutine run_case(error, bond_capacity, icase)
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in) :: bond_capacity
   integer, intent(in) :: icase
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3])
   character(len=30) :: label

   select case(icase)
   case(1, 2)
      call get_structure(mol, "MB16-43", "01")
   case(3)
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
         & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp)
   case(4)
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp)
   case(5, 6)
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
         & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp, &
         & lattice=lattice, periodic=[.true., .true., .true.])
   case(7)
      call get_structure(mol, "X23", "urea")
   case(8)
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp, &
         & lattice=lattice, periodic=[.true., .true., .true.])
   case(9)
      ! The half-cell displacement has two equivalent WS images.
      call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 3.0_wp, 0.2_wp, 0.3_wp], [3, 2]), lattice=reshape([ &
         & 6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 7.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 8.0_wp], [3, 3]), periodic=[.true., .true., .true.])
   case(10)
      ! 1D periodic wire (FTF), charged, zero nonperiodic lattice vectors
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 0.3_wp, 2.4_wp, -0.2_wp, &
         & -0.4_wp, 1.1_wp, 1.5_wp], [3, 3]), charge=1.0_wp)
      mol%lattice = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.4_wp, 5.3_wp, -0.3_wp, &
         & 0.0_wp, 0.0_wp, 0.0_wp], [3, 3])
      mol%periodic = [.false., .true., .false.]
   case(11)
      ! 2D periodic slab (TFT), charged, zero nonperiodic lattice vector
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, -0.3_wp, 0.0_wp, 2.5_wp, 0.4_wp, 0.2_wp, &
         & 0.6_wp, 0.5_wp, 2.1_wp], [3, 3]), charge=-1.0_wp)
      mol%lattice = reshape([5.2_wp, 0.1_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.5_wp, -0.2_wp, 6.1_wp], [3, 3])
      mol%periodic = [.true., .false., .true.]
   case(12)
      ! 1D periodic wire (TFF) with half-cell displacement (competing WSC images)
      call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 3.0_wp, 0.2_wp, 0.3_wp], [3, 2]))
      mol%lattice = reshape([ &
         & 6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 0.0_wp], [3, 3])
      mol%periodic = [.true., .false., .false.]
   case(13)
      ! 2D periodic slab (TFT), neutral, zero nonperiodic lattice vector
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, -0.3_wp, 0.0_wp, 2.5_wp, 0.4_wp, 0.2_wp, &
         & 0.6_wp, 0.5_wp, 2.1_wp], [3, 3]), charge=0.0_wp)
      mol%lattice = reshape([5.2_wp, 0.1_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.5_wp, -0.2_wp, 6.1_wp], [3, 3])
      mol%periodic = [.true., .false., .true.]
   case(14)
      ! 2D periodic slab (TTF) with half-cell displacement (competing WSC images)
      call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, -0.2_wp, &
         & 3.0_wp, 0.2_wp, 0.3_wp], [3, 2]), charge=0.0_wp)
      mol%lattice = reshape([ &
         & 6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 7.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 0.0_wp], [3, 3])
      mol%periodic = [.true., .true., .false.]
   case(15)
      ! 1D periodic wire (FTF), neutral, with both on-axis and off-axis pairs
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 0.16_wp, 2.12_wp, -0.12_wp, &
         & -0.4_wp, 1.1_wp, 1.5_wp], [3, 3]), charge=0.0_wp)
      mol%lattice = reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.4_wp, 5.3_wp, -0.3_wp, &
         & 0.0_wp, 0.0_wp, 0.0_wp], [3, 3])
      mol%periodic = [.false., .true., .false.]
   case default
      call test_failed(error, "Unknown charge Hessian test case")
      return
   end select
   if (bond_capacity) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call new_eeq2019_model(mol, model, error)
   end if
   if (allocated(error)) return
   ! Force an appreciable smooth-cap curvature in both models.
   if (icase == 2 .or. icase == 6) model%ncoord%cut = 2.5_wp
   call check_response(error, mol, model)
   if (allocated(error)) then
      write(label, "(a, i0, a)") "Case ", icase, ": "
      error%message = trim(label)//" "//error%message
   end if
end subroutine run_case

!> Independent reference: central differences of the weighted first response.
subroutine check_response(error, mol, model)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   class(mchrg_model_type), intent(in) :: model
   type(cache_container) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: xyz, delta, scale
   real(wp), allocatable :: q(:), weights(:), hess(:, :), numeric(:, :), again(:, :)
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :), grad(:, :), plus(:, :)
   integer :: n, i, j, a, col
   character(len=100) :: message

   n = mol%nat
   allocate(q(n), weights(n), hess(3*n, 3*n), numeric(3*n, 3*n), again(3*n, 3*n))
   allocate(dqdr(3, n, n), dqdL(3, 3, n), grad(3, n), plus(3, n))
   weights = [(sin(real(2*i, wp)), i=1, n)]
   call get_charges(model, mol, error, q, cache=cache)
   if (allocated(error)) return
   hess = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess, cache)
   if (allocated(error)) return
   do j = 1, n
      do a = 1, 3
         col = 3*(j-1)+a
         xyz = mol%xyz(a, j)
         mol%xyz(a, j) = xyz+step
         call get_charges(model, mol, error, q, dqdr, dqdL)
         if (allocated(error)) return
         plus = 0.0_wp
         do i = 1, n
            plus = plus+weights(i)*dqdr(:, :, i)
         end do
         mol%xyz(a, j) = xyz-step
         call get_charges(model, mol, error, q, dqdr, dqdL)
         if (allocated(error)) return
         grad = 0.0_wp
         do i = 1, n
            grad = grad+weights(i)*dqdr(:, :, i)
         end do
         numeric(:, col) = reshape((plus-grad)/(2*step), [3*n])
         mol%xyz(a, j) = xyz
      end do
   end do
   delta = maxval(abs(hess-numeric))
   if (delta > 1.0e-6_wp) then
      write(message, "(a, es12.4)") "Weighted charge Hessian differs from finite differences: ", delta
      call test_failed(error, trim(message))
      return
   end if
   scale = max(1.0_wp, maxval(abs(hess)))
   if (maxval(abs(hess-transpose(hess))) > 1.0e-11_wp*scale) then
      call test_failed(error, "Charge Hessian is not symmetric")
      return
   end if
   do a = 1, 3
      if (maxval(abs(sum(hess(:, a:3*n:3), dim=2))) > 1.0e-10_wp*scale) then
         call test_failed(error, "Charge Hessian violates translation invariance")
         return
      end if
   end do
   ! The same response cache may serve several weight vectors and accumulations.
   do j = 1, 3*n
      do i = 1, 3*n
         numeric(i, j) = 0.125_wp + real(i+2*j, wp)*1.0e-6_wp
      end do
   end do
   again = numeric
   call add_charge_hessian(model, mol, error, weights, again, cache)
   if (allocated(error)) return
   if (maxval(abs(again-hess-numeric)) > 1.0e-11_wp*scale) then
      call test_failed(error, "Charge Hessian does not accumulate")
      return
   end if
   again = 0.0_wp
   weights = 1.0_wp
   call add_charge_hessian(model, mol, error, weights, again, cache)
   if (allocated(error)) return
   if (maxval(abs(again)) > 1.0e-10_wp*scale) then
      call test_failed(error, "Second response violates total charge conservation")
   end if
end subroutine check_response

!> Refresh an existing cache for a different model and atom count, and exercise
!> the cache + Cartesian-response interface without a strain output.
subroutine test_cache(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(cache_container) :: cache, fresh
   class(mchrg_model_type), allocatable :: model
   real(wp) :: q1(1), q(3), refq(3), dqdr(3, 3, 3), refdr(3, 3, 3), dqdL(3, 3, 3)
   real(wp) :: hess(9, 9), refhess(9, 9), weights(3)

   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call get_charges(model, mol, error, q1, cache=cache)
   if (allocated(error)) return
   call new(mol, [30, 8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call get_charges(model, mol, error, q, dqdr=dqdr, cache=cache)
   if (allocated(error)) return
   call get_charges(model, mol, error, refq, refdr, dqdL, fresh)
   if (allocated(error)) return
   if (maxval(abs(q-refq)) > 1.0e-12_wp .or. maxval(abs(dqdr-refdr)) > 1.0e-12_wp) then
      call test_failed(error, "Refreshing the response cache changed charges or first response")
      return
   end if
   weights = [0.3_wp, -0.7_wp, 0.5_wp]
   hess = 0.0_wp
   refhess = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess, cache)
   if (allocated(error)) return
   call add_charge_hessian(model, mol, error, weights, refhess, fresh)
   if (allocated(error)) return
   if (maxval(abs(hess-refhess)) > 1.0e-12_wp) then
     call test_failed(error, "Refreshing the response cache changed second response")
   end if
end subroutine test_cache

subroutine test_invalid(error)
   type(error_type), allocatable, intent(out) :: error
   type(error_type), allocatable :: expected
   type(structure_type) :: mol
   type(cache_container) :: cache
   class(mchrg_model_type), allocatable :: model
   real(wp) :: hess(3, 3), q(1)

   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   hess = 0.125_wp
   call add_charge_hessian(model, mol, expected, [1.0_wp], hess, cache)
   if (.not. allocated(expected)) then
      call test_failed(error, "An unprepared response cache must be rejected")
      return
   end if
   call get_charges(model, mol, error, q, cache=cache)
   if (allocated(error)) return
   mol%periodic = [.true., .false., .false.]
   mol%charge = 1.0_wp
   call add_charge_hessian(model, mol, expected, [1.0_wp], hess, cache)
   if (.not. allocated(expected)) then
      call test_failed(error, "Charged 1D EEQ cell must be rejected by Hessian")
      return
   end if
   if (any(hess /= 0.125_wp)) call test_failed(error, "Rejected calls changed the Hessian")
end subroutine test_invalid

subroutine test_small_r_kernel(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: gam = 0.65_wp, alp = 0.25_wp
   real(wp) :: val_g, dr_g, drr_g, val_a, dr_a, drr_a
   real(wp) :: val_m, dr_m, drr_m, val_p, dr_p, drr_p
   real(wp) :: box_0(3, 3), box_small(3, 3), iso_lim
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp) :: cn(2), qloc(2), amat_hi(3, 3), amat_lo(3, 3), amat_0(3, 3)
   real(wp) :: dadr_hi(3, 2, 3), dadL_hi(3, 3, 3), atrace_hi(3, 2)
   real(wp) :: dadr_lo(3, 2, 3), dadL_lo(3, 3, 3), atrace_lo(3, 2)
   real(wp) :: dadr_0(3, 2, 3), dadL_0(3, 3, 3), atrace_0(3, 2)
   real(wp) :: hess_hi(6, 6), hess_lo(6, 6), hess_0(6, 6)
   real(wp), allocatable :: trans(:, :)
   real(wp), parameter :: qvec(3) = [0.35_wp, -0.35_wp, 0.0_wp]
   real(wp), parameter :: adj(2) = [0.6_wp, -0.4_wp]
   logical :: masks(3, 4)
   integer :: imask, ic

   ! 1. Direct kernel continuity at the Taylor switch (u = 1e-2) and at r = 0
   call coulomb_kernel((1.0e-2_wp - 1.0e-7_wp) / gam, gam, val_m, dr_m, drr_m)
   call coulomb_kernel((1.0e-2_wp + 1.0e-7_wp) / gam, gam, val_p, dr_p, drr_p)
   if (abs(val_m - val_p) > 1.0e-6_wp .or. abs(dr_m - dr_p) > 1.0e-6_wp &
      & .or. abs(drr_m - drr_p) > 1.0e-6_wp) then
      call test_failed(error, "Discontinuity at coulomb_kernel small-r switch")
      return
   end if

   call coulomb_kernel(0.0_wp, gam, val_g, dr_g, drr_g)
   call coulomb_kernel(0.0_wp, alp, val_a, dr_a, drr_a)
   if (abs((val_g - val_a) - 2.0_wp * (gam - alp) / sqrt(pi)) > 1.0e-15_wp) then
      call test_failed(error, "Wrong r -> 0 potential limit for Gaussian difference kernel")
      return
   end if
   if (abs(dr_g - dr_a) > 1.0e-15_wp) then
      call test_failed(error, "Wrong r -> 0 gradient limit for Gaussian difference kernel")
      return
   end if
   iso_lim = -4.0_wp * (gam**3 - alp**3) / (3.0_wp * sqrt(pi))
   box_0 = radial_hessian([0.0_wp, 0.0_wp, 0.0_wp], dr_g - dr_a, drr_g - drr_a)
   call coulomb_kernel(1.0e-8_wp, gam, val_g, dr_g, drr_g)
   call coulomb_kernel(1.0e-8_wp, alp, val_a, dr_a, drr_a)
   box_small = radial_hessian([1.0e-8_wp, 0.0_wp, 0.0_wp], dr_g - dr_a, drr_g - drr_a)
   do ic = 1, 3
      if (abs(box_0(ic, ic) - iso_lim) > 1.0e-15_wp &
         & .or. abs(box_small(ic, ic) - iso_lim) > 1.0e-13_wp) then
         call test_failed(error, "Wrong r -> 0 Hessian limit for Gaussian difference kernel")
         return
      end if
   end do

   ! 2. Continuity of EEQ potential, gradient, and Hessian across r = eps (~1.49e-8 Bohr)
   !    and down to r = 0 in 0D, 1D, 2D, and 3D
   masks = reshape([ &
      & .false., .false., .false., &
      & .true.,  .false., .false., &
      & .true.,  .true.,  .false., &
      & .true.,  .true.,  .true.], [3, 4])

   do imask = 1, 4
      ! Above eps: r = 3.0e-8 Bohr
      call new(mol, [8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, &
         & 3.0e-8_wp, 0.0_wp, 0.0_wp], [3, 2]), charge=0.0_wp, &
         & lattice=reshape([ &
         & 8.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 8.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 8.0_wp], [3, 3]), &
         & periodic=masks(:, imask))
      if (allocated(model)) deallocate(model)
      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return
      call get_dir_trans(mol, trans, 25.0_wp)
      call model%ncoord%get_coordination_number(mol, trans, cn)
      call model%local_charge(mol, trans, qloc)
      call model%update(mol, cache, cn, qloc)
      allocate(cache%qvec(3), source=qvec)
      select type (raw => cache%raw)
      type is (eeq_cache)
         allocate(raw%dcndr(3, 2, 2), source=0.0_wp)
      end select
      call model%get_coulomb_matrix(mol, cache, amat_hi)
      call model%get_coulomb_derivs(mol, cache, qvec, dadr_hi, dadL_hi, atrace_hi)
      hess_hi = 0.0_wp
      call model%add_model_hessian(mol, error, cache, adj, hess_hi)
      if (allocated(error)) return
      deallocate(cache%qvec, cache%raw)

      ! Below eps: r = 5.0e-9 Bohr (previously skipped by r < eps!)
      mol%xyz(:, 2) = [5.0e-9_wp, 0.0_wp, 0.0_wp]
      call model%update(mol, cache, cn, qloc)
      allocate(cache%qvec(3), source=qvec)
      select type (raw => cache%raw)
      type is (eeq_cache)
         allocate(raw%dcndr(3, 2, 2), source=0.0_wp)
      end select
      call model%get_coulomb_matrix(mol, cache, amat_lo)
      call model%get_coulomb_derivs(mol, cache, qvec, dadr_lo, dadL_lo, atrace_lo)
      hess_lo = 0.0_wp
      call model%add_model_hessian(mol, error, cache, adj, hess_lo)
      if (allocated(error)) return
      deallocate(cache%qvec, cache%raw)

      ! At r = 1.0e-12 Bohr -> 0
      mol%xyz(:, 2) = [1.0e-12_wp, 0.0_wp, 0.0_wp]
      call model%update(mol, cache, cn, qloc)
      allocate(cache%qvec(3), source=qvec)
      select type (raw => cache%raw)
      type is (eeq_cache)
         allocate(raw%dcndr(3, 2, 2), source=0.0_wp)
      end select
      call model%get_coulomb_matrix(mol, cache, amat_0)
      call model%get_coulomb_derivs(mol, cache, qvec, dadr_0, dadL_0, atrace_0)
      hess_0 = 0.0_wp
      call model%add_model_hessian(mol, error, cache, adj, hess_0)
      if (allocated(error)) return
      deallocate(cache%qvec, cache%raw)

      if (maxval(abs(amat_hi - amat_lo)) > 1.0e-12_wp &
         & .or. maxval(abs(amat_lo - amat_0)) > 1.0e-12_wp) then
         call test_failed(error, "Coulomb potential discontinuous as off-diagonal r -> 0")
         return
      end if
      if (maxval(abs(dadr_hi - dadr_lo)) > 1.0e-7_wp &
         & .or. maxval(abs(dadr_0)) > 1.0e-10_wp) then
         call test_failed(error, "Coulomb gradient discontinuous as off-diagonal r -> 0")
         return
      end if
      if (maxval(abs(hess_hi - hess_lo)) > 1.0e-12_wp &
         & .or. maxval(abs(hess_lo - hess_0)) > 1.0e-12_wp) then
         call test_failed(error, "Coulomb Hessian discontinuous as off-diagonal r -> 0")
         return
      end if
   end do

end subroutine test_small_r_kernel

end module test_hessian
