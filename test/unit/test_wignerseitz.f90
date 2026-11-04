! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module test_wignerseitz
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io_structure, only: structure_type, new
   use mctc_wignerseitz, only: wignerseitz_cell, new_wignerseitz_cell, &
      & get_wignerseitz_weights
   use mstore, only: get_structure
   use multicharge, only: mchrg_model_type, new_eeq2019_model, new_eeqbc2025_model, &
      & get_charges
   implicit none
   private

   public :: collect_wignerseitz

   real(wp), parameter :: thr = 100 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_wignerseitz(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("lattice-points-0d", test_latticepoints_0d), &
      & new_unittest("lattice-points-1d", test_latticepoints_1d), &
      & new_unittest("lattice-points-2d", test_latticepoints_2d), &
      & new_unittest("lattice-points-3d", test_latticepoints_3d), &
      & new_unittest("wignerseitz-3d-skewed-counterexample", test_wsc_3d_skewed_counterexample), &
      & new_unittest("wignerseitz-65deg-counterexample", test_wsc_65deg_counterexample), &
      & new_unittest("wignerseitz-cell-0d", test_wsc_cell_0d), &
      & new_unittest("wignerseitz-cell-1d", test_wsc_cell_1d), &
      & new_unittest("wignerseitz-cell-2d", test_wsc_cell_2d), &
      & new_unittest("wignerseitz-cell-3d", test_wsc_cell_3d), &
      & new_unittest("wignerseitz-resource-limits", test_wsc_resource_limits), &
      & new_unittest("wignerseitz-robustness-1d", test_wsc_robustness_1d), &
      & new_unittest("wignerseitz-robustness-2d", test_wsc_robustness_2d), &
      & new_unittest("wignerseitz-robustness-2d-corner", test_wsc_robustness_2d_corner), &
      & new_unittest("wignerseitz-self-images", test_wsc_self_images) &
      & ]

end subroutine collect_wignerseitz


subroutine test_latticepoints_0d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: trans(:, :)

   real(wp), parameter :: cutoff = 25.0_wp

   call get_structure(mol, "MB16-43", "01")

   call get_lattice_points(mol%periodic, mol%lattice, thr2, trans, error)
   if (allocated(error)) return

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 1)
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 1)
   if (allocated(error)) return

end subroutine test_latticepoints_0d


subroutine test_latticepoints_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: trans1(:, :), trans2(:, :)
   real(wp) :: lat(3, 3)
   integer :: imask, itr

   integer, parameter :: nat = 1
   integer, parameter :: num(nat) = [6]
   real(wp), parameter :: a1(3) = [5.0_wp, 0.0_wp, 0.0_wp]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, nat])
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .false., .false., &
      & .false., .true., .false., &
      & .false., .false., .true.], &
      & [3, 3])

   do imask = 1, 3
      lat = 0.0_wp
      lat(:, imask) = a1
      call new(mol, num, xyz, lattice=lat, periodic=masks(:, imask))

      call get_lattice_points(mol%periodic, mol%lattice, thr2, trans1, error)
      if (allocated(error)) return

      call check(error, size(trans1, 1), 3)
      if (allocated(error)) return

      call check(error, size(trans1, 2), 3)
      if (allocated(error)) return

      call get_lattice_points(mol%periodic, mol%lattice, 12.0_wp, trans1, error)
      if (allocated(error)) return

      ! ceiling(12.0 / 5.0) = 3 -> 2*3 + 1 = 7 points
      call check(error, size(trans1, 2), 7)
      if (allocated(error)) return

      ! Perturb nonperiodic lattice vectors wildly; translations must be identical
      do itr = 1, 3
         if (itr /= imask) then
            lat(:, itr) = [100.0_wp * itr, -50.0_wp, 250.0_wp]
         end if
      end do
      mol%lattice = lat
      call get_lattice_points(mol%periodic, mol%lattice, 12.0_wp, trans2, error)
      if (allocated(error)) return

      call check(error, size(trans2, 2), 7)
      if (allocated(error)) return

      if (maxval(abs(trans1 - trans2)) > thr) then
         call test_failed(error, "1D lattice points depend on nonperiodic lattice vectors")
         return
      end if
   end do

end subroutine test_latticepoints_1d


subroutine test_latticepoints_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: trans1(:, :), trans2(:, :)
   real(wp) :: lat(3, 3)
   integer :: imask, p1, p2, np

   integer, parameter :: nat = 1
   integer, parameter :: num(nat) = [6]
   real(wp), parameter :: a1(3) = [5.0_wp, 0.0_wp, 0.0_wp]
   real(wp), parameter :: a2(3) = [0.0_wp, 6.0_wp, 0.0_wp]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, nat])
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .true., .false., &
      & .true., .false., .true., &
      & .false., .true., .true.], &
      & [3, 3])

   do imask = 1, 3
      select case (imask)
      case (1)
         p1 = 1
         p2 = 2
         np = 3
      case (2)
         p1 = 1
         p2 = 3
         np = 2
      case default
         p1 = 2
         p2 = 3
         np = 1
      end select

      lat = 0.0_wp
      lat(:, p1) = a1
      lat(:, p2) = a2
      call new(mol, num, xyz, lattice=lat, periodic=masks(:, imask))

      call get_lattice_points(mol%periodic, mol%lattice, thr2, trans1, error)
      if (allocated(error)) return

      call check(error, size(trans1, 2), 9)
      if (allocated(error)) return

      call get_lattice_points(mol%periodic, mol%lattice, 11.0_wp, trans1, error)
      if (allocated(error)) return

      ! rep1 = ceiling(11/5) = 3 (7), rep2 = ceiling(11/6) = 2 (5) -> 35
      call check(error, size(trans1, 2), 35)
      if (allocated(error)) return

      lat(:, np) = [42.0_wp, -99.0_wp, 500.0_wp]
      mol%lattice = lat
      call get_lattice_points(mol%periodic, mol%lattice, 11.0_wp, trans2, error)
      if (allocated(error)) return

      if (maxval(abs(trans1 - trans2)) > thr) then
         call test_failed(error, "2D lattice points depend on nonperiodic lattice vector")
         return
      end if
   end do

end subroutine test_latticepoints_2d


subroutine test_latticepoints_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: trans(:, :)

   real(wp), parameter :: cutoff = 25.0_wp

   call get_structure(mol, "ICE10", "ii")

   call get_lattice_points(mol%periodic, mol%lattice, thr2, trans, error)
   if (allocated(error)) return

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 27)
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 343)
   if (allocated(error)) return

end subroutine test_latticepoints_3d


subroutine test_wsc_3d_skewed_counterexample(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc
   real(wp) :: rij(3), vec(3), d2_found, weights(16)
   integer :: img

   real(wp), parameter :: d2_expect = 6.5_wp
   integer, parameter :: nat = 2
   integer, parameter :: num(nat) = [6, 8]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      &  0.0_wp, 0.0_wp, 0.0_wp, &
      & 12.5_wp, 0.5_wp, 0.0_wp], &
      & [3, nat])
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      &  5.0_wp, 0.0_wp, 0.0_wp, &
      & 25.0_wp, 1.0_wp, 0.0_wp, &
      &  0.0_wp, 0.0_wp, 5.0_wp], &
      & [3, 3])
   logical, parameter :: periodic(3) = [.true., .true., .true.]

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   call new_wignerseitz_cell(wsc, mol, error)
   if (allocated(error)) return

   call check(error, wsc%nimg(1, 2), 4)
   if (allocated(error)) return

   call check(error, wsc%nimg(2, 1), 4)
   if (allocated(error)) return

   rij = mol%xyz(:, 2) - mol%xyz(:, 1)
   do img = 1, wsc%nimg(1, 2)
      vec = rij - wsc%trans(:, wsc%tridx(img, 1, 2))
      d2_found = sum(vec**2)
      if (abs(d2_found - d2_expect) > thr) then
         call test_failed(error, "3D skewed counterexample did not find d^2 = 6.5 bohr^2")
         return
      end if
   end do

   call get_wignerseitz_weights(wsc, 1, 2, rij, weights)
   if (maxval(abs(weights(1:4) - 0.25_wp)) > thr) then
      call test_failed(error, "3D skewed counterexample weights are not [0.25, 0.25, 0.25, 0.25]")
      return
   end if

end subroutine test_wsc_3d_skewed_counterexample


subroutine test_wsc_65deg_counterexample(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc
   real(wp) :: rij(3), vec1(3), vec2(3)
   real(wp) :: d2_found1, d2_found2, d2_expect, weights(16)

   real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: angle = 65.0_wp * pi / 180.0_wp
   integer, parameter :: nat = 2
   integer, parameter :: num(nat) = [6, 8]
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 5.0_wp, 0.0_wp, 0.0_wp, &
      & 5.0_wp * cos(angle), 5.0_wp * sin(angle), 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, 3])
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.39_wp * (lattice(:, 1) + lattice(:, 2))], &
      & [3, nat])
   logical, parameter :: periodic(3) = [.true., .true., .false.]

   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   call new_wignerseitz_cell(wsc, mol, error)
   if (allocated(error)) return

   call check(error, wsc%nimg(1, 2), 2)
   if (allocated(error)) return

   call check(error, wsc%nimg(2, 1), 2)
   if (allocated(error)) return

   rij = mol%xyz(:, 2) - mol%xyz(:, 1)
   d2_expect = sum((rij - 2.0_wp * lattice(:, 1) - lattice(:, 2))**2)
   vec1 = rij - wsc%trans(:, wsc%tridx(1, 1, 2))
   vec2 = rij - wsc%trans(:, wsc%tridx(2, 1, 2))
   d2_found1 = sum(vec1**2)
   d2_found2 = sum(vec2**2)

   if (abs(d2_found1 - d2_expect) > thr .or. abs(d2_found2 - d2_expect) > thr) then
      call test_failed(error, "65-degree 2D counterexample failed to find (2,1) and (1,2) images")
      return
   end if

   call get_wignerseitz_weights(wsc, 1, 2, rij, weights)
   if (maxval(abs(weights(1:2) - 0.5_wp)) > thr) then
      call test_failed(error, "65-degree 2D counterexample weights are not [0.5, 0.5]")
      return
   end if

end subroutine test_wsc_65deg_counterexample


subroutine test_wsc_cell_0d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc

   call get_structure(mol, "MB16-43", "02")

   call new_wignerseitz_cell(wsc, mol, error)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 2), 1)
   if (allocated(error)) return

end subroutine test_wsc_cell_0d


subroutine test_wsc_cell_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc
   real(wp) :: lat(3, 3), vec(3), total(3)
   integer :: imask, i

   integer, parameter :: nat = 2
   integer, parameter :: num(nat) = [8, 1]
   real(wp), parameter :: a1(3) = [3.0_wp, 4.0_wp, 0.0_wp]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.5_wp, 2.0_wp, 1.2_wp], &
      & [3, nat])
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .false., .false., &
      & .false., .true., .false., &
      & .false., .false., .true.], &
      & [3, 3])

   do imask = 1, 3
      lat = 0.0_wp
      lat(:, imask) = a1
      call new(mol, num, xyz, lattice=lat, periodic=masks(:, imask))
      call new_wignerseitz_cell(wsc, mol, error)
      if (allocated(error)) return

      call check(error, size(wsc%trans, 2), 3)
      if (allocated(error)) return

      ! Atom 1 has 2 closest self images (+a and -a)
      call check(error, wsc%nimg(1, 1), 2)
      if (allocated(error)) return

      ! Atom 2 is at half-cell along a, so it has 2 equidistant images to atom 1
      call check(error, wsc%nimg(2, 1), 2)
      if (allocated(error)) return

      total = 0.0_wp
      do i = 1, wsc%nimg(1, 1)
         vec = wsc%trans(:, wsc%tridx(i, 1, 1))
         if (abs(norm2(vec) - 5.0_wp) > thr) then
            call test_failed(error, "1D WSC self image has wrong length")
            return
         end if

         total = total + vec
      end do

      if (maxval(abs(total)) > thr) then
         call test_failed(error, "1D WSC self images do not cancel")
         return
      end if
   end do

end subroutine test_wsc_cell_1d


subroutine test_wsc_cell_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc
   real(wp) :: total(3), vec(3)
   integer :: i

   integer, parameter :: nat = 1
   integer, parameter :: num(nat) = [8]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, nat])
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 1.5_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 5.2_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 5.2_wp], &
      & [3, 3])
   logical, parameter :: periodic(3) = [.false., .true., .true.]

   ! Square 2D lattice in YZ plane (FTT) with nonperiodic X vector shorter than Y/Z!
   ! A broken 3D WSC would pick the X vector as closest self image; 2D WSC must ignore X.
   call new(mol, num, xyz, lattice=lattice, periodic=periodic)
   call new_wignerseitz_cell(wsc, mol, error)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 2), 9)
   if (allocated(error)) return

   call check(error, wsc%nimg(1, 1), 4)
   if (allocated(error)) return

   total = 0.0_wp
   do i = 1, wsc%nimg(1, 1)
      vec = wsc%trans(:, wsc%tridx(i, 1, 1))
      if (abs(vec(1)) > thr .or. abs(norm2(vec) - 5.2_wp) > thr) then
         call test_failed(error, "2D WSC selected nonperiodic or non-minimal self image")
         return
      end if

      total = total + vec
   end do

   if (maxval(abs(total)) > thr) then
      call test_failed(error, "2D WSC self images do not cancel")
      return
   end if

end subroutine test_wsc_cell_2d


subroutine test_wsc_cell_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc

   call get_structure(mol, "X23", "pyrazole")

   call new_wignerseitz_cell(wsc, mol, error)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 2), 105)
   if (allocated(error)) return

end subroutine test_wsc_cell_3d


subroutine test_wsc_resource_limits(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc
   type(error_type), allocatable :: call_err
   class(mchrg_model_type), allocatable :: model_eeq, model_eeqbc
   real(wp) :: qvec(2)
   real(wp), allocatable :: trans(:, :)

   integer, parameter :: nat = 2
   integer, parameter :: num(nat) = [6, 8]
   real(wp), parameter :: lat_1d(3, 3) = reshape([ &
      & 5.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, 3])
   real(wp), parameter :: xyz_1d(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 500001.2_wp, 0.0_wp, 0.0_wp], &
      & [3, nat])
   real(wp), parameter :: lat_2d(3, 3) = reshape([ &
      & 5.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 5.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, 3])
   real(wp), parameter :: xyz_2d(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1501.2_wp, 1501.3_wp, 0.0_wp], &
      & [3, nat])
   real(wp), parameter :: lat_3d(3, 3) = reshape([ &
      & 5.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 5.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 5.0_wp], &
      & [3, 3])
   logical, parameter :: pbc_1d(3) = [.true., .false., .false.]
   logical, parameter :: pbc_2d(3) = [.true., .true., .false.]

   ! 1. 1D resource-limit regression (> 100,000 candidate periods)
   call new(mol, num, xyz_1d, lattice=lat_1d, periodic=pbc_1d)

   call new_wignerseitz_cell(wsc, mol, call_err)
   if (.not. allocated(call_err) .or. allocated(wsc%trans)) then
      call test_failed(error, "1D resource limit did not fail cleanly in new_wignerseitz_cell")
      return
   end if

   call new_eeq2019_model(mol, model_eeq, error)
   if (allocated(error)) return

   call get_charges(model_eeq, mol, call_err, qvec)
   if (.not. allocated(call_err)) then
      call test_failed(error, "1D WSC resource limit was not propagated by EEQ get_charges")
      return
   end if

   call new_eeqbc2025_model(mol, model_eeqbc, error)
   if (allocated(error)) return

   call get_charges(model_eeqbc, mol, call_err, qvec)
   if (.not. allocated(call_err)) then
      call test_failed(error, "1D WSC resource limit was not propagated by EEQBC get_charges")
      return
   end if

   ! 2. 2D resource-limit regression (> 100,000 candidate translations)
   call new(mol, num, xyz_2d, lattice=lat_2d, periodic=pbc_2d)

   call new_wignerseitz_cell(wsc, mol, call_err)
   if (.not. allocated(call_err) .or. allocated(wsc%trans)) then
      call test_failed(error, "2D resource limit did not fail cleanly in new_wignerseitz_cell")
      return
   end if

   call get_charges(model_eeq, mol, call_err, qvec)
   if (.not. allocated(call_err)) then
      call test_failed(error, "2D WSC resource limit was not propagated by EEQ get_charges")
      return
   end if

   call get_charges(model_eeqbc, mol, call_err, qvec)
   if (.not. allocated(call_err)) then
      call test_failed(error, "2D WSC resource limit was not propagated by EEQBC get_charges")
      return
   end if

   ! 3. Lattice point enumeration overflow (rep = [650, 650, 650]) and resource limit
   call get_lattice_points(lat_3d, [650, 650, 650], .true., trans, call_err)
   if (.not. allocated(call_err) .or. allocated(trans)) then
      call test_failed(error, "rep = [650, 650, 650] did not fail cleanly in get_lattice_points")
      return
   end if

   call get_lattice_points(lat_3d, [200, 200, 200], .true., trans, call_err)
   if (.not. allocated(call_err) .or. allocated(trans)) then
      call test_failed(error, "rep = [200, 200, 200] did not fail cleanly in get_lattice_points")
      return
   end if

end subroutine test_wsc_resource_limits


subroutine test_wsc_robustness_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol_ref, mol_test
   type(wignerseitz_cell) :: wsc_ref, wsc_test
   real(wp) :: v_ref(3), v_test(3)
   integer :: iat, jat

   integer, parameter :: nat = 3
   integer, parameter :: num(nat) = [6, 8, 1]
   real(wp), parameter :: a1(3) = [3.0_wp, 4.0_wp, 12.0_wp]
   real(wp), parameter :: xyz_ref(3, nat) = reshape([ &
      &  0.1_wp,  0.2_wp, -0.3_wp, &
      &  1.1_wp,  1.4_wp,  3.5_wp, &
      & -0.8_wp, -0.9_wp, -2.5_wp], &
      & [3, nat])
   real(wp), parameter :: xyz_test(3, nat) = reshape([ &
      & xyz_ref(:, 1), &
      & xyz_ref(:, 2) + 5.0_wp * a1, &
      & xyz_ref(:, 3) - 4.0_wp * a1], &
      & [3, nat])
   real(wp), parameter :: lat_ref(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 3.0_wp, 4.0_wp, 12.0_wp, &
      & 0.0_wp, 0.0_wp,  0.0_wp], &
      & [3, 3])
   logical, parameter :: periodic(3) = [.false., .true., .false.]

   ! 1D wire rotated in 3D space with zero nonperiodic lattice columns
   ! and multi-cell atom translations (|n| > 1)
   call new(mol_ref, num, xyz_ref, lattice=lat_ref, periodic=periodic)
   call new(mol_test, num, xyz_test, lattice=lat_ref, periodic=periodic)
   call new_wignerseitz_cell(wsc_ref, mol_ref, error)
   if (allocated(error)) return

   call new_wignerseitz_cell(wsc_test, mol_test, error)
   if (allocated(error)) return

   do iat = 1, 3
      do jat = 1, 3
         call check(error, wsc_test%nimg(jat, iat), wsc_ref%nimg(jat, iat))
         if (allocated(error)) return

         if (iat /= jat) then
            v_ref = mol_ref%xyz(:, jat) - mol_ref%xyz(:, iat) &
               & + wsc_ref%trans(:, wsc_ref%tridx(1, jat, iat))
            v_test = mol_test%xyz(:, jat) - mol_test%xyz(:, iat) &
               & + wsc_test%trans(:, wsc_test%tridx(1, jat, iat))
            if (maxval(abs(v_test - v_ref)) > thr) then
               call test_failed(error, "1D WSC failed for multi-cell atom displacement")
               return
            end if
         end if
      end do
   end do

end subroutine test_wsc_robustness_1d


subroutine test_wsc_robustness_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol_ref, mol_test
   type(wignerseitz_cell) :: wsc_ref, wsc_test
   real(wp) :: v_ref(3), v_test(3)
   integer :: iat, jat

   integer, parameter :: nat = 3
   integer, parameter :: num(nat) = [6, 8, 1]
   real(wp), parameter :: a1(3) = [4.0_wp, 3.0_wp, 0.0_wp]
   real(wp), parameter :: a2(3) = [-3.0_wp, 4.0_wp, 5.0_wp]
   real(wp), parameter :: xyz_ref(3, nat) = reshape([ &
      &  0.2_wp, -0.1_wp, 0.3_wp, &
      &  1.4_wp,  1.2_wp, 1.1_wp, &
      & -1.1_wp,  0.8_wp, 1.7_wp], &
      & [3, nat])
   real(wp), parameter :: xyz_test(3, nat) = reshape([ &
      & xyz_ref(:, 1), &
      & xyz_ref(:, 2) + 3.0_wp * a1 - 4.0_wp * a2, &
      & xyz_ref(:, 3) - 5.0_wp * a1 + 2.0_wp * a2], &
      & [3, nat])
   real(wp), parameter :: lat_ref(3, 3) = reshape([ &
      &  4.0_wp, 3.0_wp, 0.0_wp, &
      &  0.0_wp, 0.0_wp, 0.0_wp, &
      & -3.0_wp, 4.0_wp, 5.0_wp], &
      & [3, 3])
   real(wp), parameter :: lat_test(3, 3) = reshape([ &
      & a1, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & a2 + 6.0_wp * a1], &
      & [3, 3])
   logical, parameter :: periodic(3) = [.true., .false., .true.]

   ! 2D slab in rotated 3D subspace with zero nonperiodic column,
   ! skewed unimodular basis transformation (a2' = a2 + 6*a1), and multi-cell shifts
   call new(mol_ref, num, xyz_ref, lattice=lat_ref, periodic=periodic)
   call new(mol_test, num, xyz_test, lattice=lat_test, periodic=periodic)
   call new_wignerseitz_cell(wsc_ref, mol_ref, error)
   if (allocated(error)) return

   call new_wignerseitz_cell(wsc_test, mol_test, error)
   if (allocated(error)) return

   do iat = 1, 3
      do jat = 1, 3
         call check(error, wsc_test%nimg(jat, iat), wsc_ref%nimg(jat, iat))
         if (allocated(error)) return

         if (iat /= jat) then
            v_ref = mol_ref%xyz(:, jat) - mol_ref%xyz(:, iat) &
               & + wsc_ref%trans(:, wsc_ref%tridx(1, jat, iat))
            v_test = mol_test%xyz(:, jat) - mol_test%xyz(:, iat) &
               & + wsc_test%trans(:, wsc_test%tridx(1, jat, iat))
            if (maxval(abs(v_test - v_ref)) > thr) then
               call test_failed(error, &
                  & "2D WSC failed under skewed unimodular basis and multi-cell shift")
               return
            end if
         end if
      end do
   end do

end subroutine test_wsc_robustness_2d


subroutine test_wsc_robustness_2d_corner(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol_ref, mol_test
   type(wignerseitz_cell) :: wsc_ref, wsc_test
   real(wp) :: v_test(3)

   integer, parameter :: nat = 2
   integer, parameter :: num(nat) = [8, 1]
   real(wp), parameter :: xyz_ref(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp, 3.0_wp, 0.4_wp], &
      & [3, nat])
   real(wp), parameter :: xyz_test(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp + 0.03_wp, 3.0_wp - 0.03_wp, 0.4_wp], &
      & [3, nat])
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 6.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 6.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, 3])
   logical, parameter :: periodic(3) = [.true., .true., .false.]

   ! Exact corner of 2D WS cell: (3.0, 3.0, 0.4) -> 4 equidistant images
   call new(mol_ref, num, xyz_ref, lattice=lattice, periodic=periodic)
   call new_wignerseitz_cell(wsc_ref, mol_ref, error)
   if (allocated(error)) return

   call check(error, wsc_ref%nimg(2, 1), 4)
   if (allocated(error)) return

   ! Perturb by +0.03 Bohr across the WS boundary (outside 0.3 Bohr^2 tie window)
   call new(mol_test, num, xyz_test, lattice=lattice, periodic=periodic)
   call new_wignerseitz_cell(wsc_test, mol_test, error)
   if (allocated(error)) return

   call check(error, wsc_test%nimg(2, 1), 1)
   if (allocated(error)) return

   v_test = mol_test%xyz(:, 2) - mol_test%xyz(:, 1) &
      & + wsc_test%trans(:, wsc_test%tridx(1, 2, 1))
   if (abs(v_test(1) - (-3.0_wp + 0.03_wp)) > thr &
      & .or. abs(v_test(2) - (3.0_wp - 0.03_wp)) > thr) then
      call test_failed(error, "2D WSC selected wrong image after boundary perturbation")
      return
   end if

end subroutine test_wsc_robustness_2d_corner


!> Excluding the zero image must preserve the original translation indices.
subroutine test_wsc_self_images(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell) :: wsc
   real(wp) :: lattice(3, 3), vec(3), total(3)
   integer :: icase, i, expected

   integer, parameter :: nat = 1
   integer, parameter :: num(nat) = [8]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp], &
      & [3, nat])
   real(wp), parameter :: lat_skew(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, &
      & 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], &
      & [3, 3])
   real(wp), parameter :: lat_cubic(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 5.1_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 5.1_wp], &
      & [3, 3])
   logical, parameter :: periodic(3) = [.true., .true., .true.]

   do icase = 1, 2
      lattice = lat_skew
      expected = 2
      if (icase == 2) then
         lattice = lat_cubic
         expected = 6
      end if

      call new(mol, num, xyz, lattice=lattice, periodic=periodic)
      call new_wignerseitz_cell(wsc, mol, error)
      if (allocated(error)) return

      call check(error, wsc%nimg(1, 1), expected)
      if (allocated(error)) return

      call check(error, wsc%nimg_max, expected)
      if (allocated(error)) return

      total = 0.0_wp
      do i = 1, wsc%nimg(1, 1)
         vec = wsc%trans(:, wsc%tridx(i, 1, 1))
         if (abs(norm2(vec) - 5.1_wp) > thr) then
            call test_failed(error, "Selected self image is not a closest lattice point")
            return
         end if

         total = total + vec
      end do

      if (maxval(abs(total)) > thr) then
         call test_failed(error, "Closest self images do not occur in opposite pairs")
         return
      end if
   end do

end subroutine test_wsc_self_images


end module test_wignerseitz
