program controlled_love_numbers

  use, intrinsic :: ieee_arithmetic, only : ieee_divide_by_zero, &
      ieee_get_flag, ieee_invalid, ieee_is_finite, ieee_set_flag
  use, intrinsic :: iso_fortran_env, only : error_unit
  use module_constants
  use module_DECK
  use module_force
  use module_matrix
  use module_physical_constants, only : acceleration_norm, &
                                        density_norm, &
                                        gravitational_potential_norm, &
                                        length_norm, load_norm
  use module_PREM
  use module_spherical_mesh
  use full_domain_matrix
  implicit none

  character(len=512) :: argument,model_file
  integer(i4b) :: argument_index,element_count,number_of_degrees
  integer(i4b), dimension(:), allocatable :: degrees
  logical :: divide_by_zero,invalid
  real(dp) :: maximum_radial_step,surface_gravity,surface_radius
  type(spherical_model) :: model
  type(spherical_model_mesh) :: mesh

  number_of_degrees = command_argument_count()-2
  if(number_of_degrees < 1) then
     error stop 'usage: controlled_love_numbers model drmax degree...'
  end if

  call get_command_argument(1,model_file)
  call get_command_argument(2,argument)
  read(argument,*) maximum_radial_step
  if(maximum_radial_step <= 0.0_dp) error stop 'drmax must be positive'

  allocate(degrees(number_of_degrees))
  do argument_index = 1,number_of_degrees
     call get_command_argument(argument_index+2,argument)
     read(argument,*) degrees(argument_index)
     if(degrees(argument_index) < 1) error stop 'degrees must be positive'
  end do

  call ieee_set_flag(ieee_invalid,.false.)
  call ieee_set_flag(ieee_divide_by_zero,.false.)

  if(trim(model_file) == 'analytic-prem') then
     model = elastic_PREM(.false.)
  else
     model = elastic_DECK(trim(model_file))
  end if
  mesh = spherical_mesh(5,model,maximum_radial_step)
  element_count = count_elements(mesh)

  associate(layer => mesh%section(mesh%nsections)%layer( &
      mesh%section(mesh%nsections)%nlayers))
    surface_radius = layer%r(layer%ngll,layer%nspec)*length_norm
    surface_gravity = layer%g(layer%ngll,layer%nspec)*acceleration_norm
  end associate
  write(error_unit,'(a,1x,2(es26.17e3,1x),2(i0,1x))') &
      'metadata',surface_radius,surface_gravity,element_count, &
      4*element_count+1
  call write_layer_diagnostics(mesh)
  call ieee_set_flag(ieee_invalid,.false.)
  call ieee_set_flag(ieee_divide_by_zero,.false.)

  do argument_index = 1,number_of_degrees
     call solve_degree(mesh,degrees(argument_index))
  end do
  call ieee_get_flag(ieee_invalid,invalid)
  call ieee_get_flag(ieee_divide_by_zero,divide_by_zero)
  write(error_unit,'(a,1x,2(l1,1x))') &
      'floating_point_flags',invalid,divide_by_zero
  if(invalid .or. divide_by_zero) then
     error stop 'controlled calculation raised invalid or divide by zero'
  end if

contains

  integer(i4b) function count_elements(mesh) result(count)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b) :: isection,ilayer

    count = 0
    do isection = 1,mesh%nsections
       do ilayer = 1,mesh%section(isection)%nlayers
          count = count+mesh%section(isection)%layer(ilayer)%nspec
       end do
    end do
  end function count_elements


  subroutine write_layer_diagnostics(mesh)
    type(spherical_model_mesh), intent(in) :: mesh

    integer(i4b) :: flat_layer,isection,ilayer,number_of_layers

    number_of_layers = 0
    do isection = 1,mesh%nsections
       number_of_layers = number_of_layers+mesh%section(isection)%nlayers
    end do
    write(error_unit,'(a,1x,i0)',advance='no') &
        'layers',number_of_layers
    do isection = 1,mesh%nsections
       do ilayer = 1,mesh%section(isection)%nlayers
          write(error_unit,'(1x,i0)',advance='no') &
              mesh%section(isection)%layer(ilayer)%nspec
       end do
    end do
    write(error_unit,*)

    flat_layer = 0
    do isection = 1,mesh%nsections
       do ilayer = 1,mesh%section(isection)%nlayers
          associate(layer => mesh%section(isection)%layer(ilayer))
            select type(layer)
            class is(spherical_fluid_elastic_layer_mesh)
              write(error_unit,'(a,1x,i0,1x,4(es26.17e3,1x))') &
                  'fluid',flat_layer, &
                  layer%rho(1,1)*density_norm, &
                  layer%rho(layer%ngll,layer%nspec)*density_norm, &
                  minval(layer%drho)*density_norm/length_norm, &
                  maxval(layer%drho)*density_norm/length_norm
            end select
          end associate
          flat_layer = flat_layer+1
       end do
    end do
  end subroutine write_layer_diagnostics


  subroutine solve_degree(mesh,degree)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: degree

    integer(i4b) :: info,surface_element,surface_layer,surface_node
    integer(i4b) :: surface_potential,surface_radial
    real(dp) :: h_phi,h_t,h_u,k_phi,k_t,k_u,reciprocity_error
    real(dp) :: surface_gravity,surface_radius
    real(dp), dimension(3) :: residuals
    real(dp), dimension(:,:), allocatable :: original_matrix
    real(dp), dimension(:,:), allocatable :: right_hand_sides
    real(dp), dimension(:,:), allocatable :: right_hand_side_copy
    type(radial_matrix) :: matrix

    matrix = build_spheroidal_matrix_from_radius( &
        mesh,degree,0.0_dp,.false.)
    if(.not.all(ieee_is_finite(matrix%a))) &
        error stop 'controlled matrix contains a non-finite entry'
    allocate(original_matrix(matrix%ldab,matrix%ndim))
    original_matrix = matrix%a
    allocate(right_hand_sides(matrix%ndim,3))
    allocate(right_hand_side_copy(matrix%ndim,3))
    right_hand_sides = 0.0_dp

    surface_layer = mesh%section(mesh%nsections)%nlayers
    associate(layer => mesh%section(mesh%nsections)%layer(surface_layer), &
              ibool => matrix%ibool%section(mesh%nsections)%layer( &
                  surface_layer))
      select type(layer)
      class is(spherical_solid_elastic_layer_mesh)
        surface_node = layer%ngll
        surface_element = layer%nspec
        surface_radius = layer%r(surface_node,surface_element)
        surface_gravity = layer%g(surface_node,surface_element)
        surface_radial = ibool%get(1,surface_node,surface_element)

        right_hand_sides(surface_radial,1) = &
            -surface_gravity*surface_radius*surface_radius
        if(degree > 1) then
           surface_potential = ibool%get(3,surface_node,surface_element)
           right_hand_sides(surface_potential,2) = &
               -surface_radius*surface_radius
        else
           surface_potential = 0
        end if
      class default
        error stop 'controlled model requires a solid surface'
      end select
    end associate

    call force_for_unit_harmonic_tide( &
        mesh,matrix%ibool,degree,right_hand_sides(:,3:3))
    if(.not.all(ieee_is_finite(right_hand_sides))) &
        error stop 'controlled right-hand side contains a non-finite entry'
    right_hand_side_copy = right_hand_sides

    call dpbtrf('U',matrix%ndim,matrix%kd,matrix%a,matrix%ldab,info)
    if(info /= 0) error stop 'controlled matrix factorisation failed'
    call dpbtrs('U',matrix%ndim,matrix%kd,3,matrix%a,matrix%ldab, &
                right_hand_sides,matrix%ndim,info)
    if(info /= 0) error stop 'controlled three-column solve failed'
    if(.not.all(ieee_is_finite(right_hand_sides))) &
        error stop 'controlled solution contains a non-finite entry'

    call relative_residual(original_matrix,right_hand_sides(:,1), &
                           right_hand_side_copy(:,1),matrix%kd,residuals(1))
    call relative_residual(original_matrix,right_hand_sides(:,2), &
                           right_hand_side_copy(:,2),matrix%kd,residuals(2))
    call relative_residual(original_matrix,right_hand_sides(:,3), &
                           right_hand_side_copy(:,3),matrix%kd,residuals(3))

    h_u = right_hand_sides(surface_radial,1)*length_norm/load_norm
    h_phi = right_hand_sides(surface_radial,2)*length_norm/load_norm
    h_t = right_hand_sides(surface_radial,3)*length_norm
    if(degree > 1) then
       k_u = right_hand_sides(surface_potential,1) * &
             gravitational_potential_norm/load_norm
       k_phi = right_hand_sides(surface_potential,2) * &
               gravitational_potential_norm/load_norm
       k_t = right_hand_sides(surface_potential,3) * &
             gravitational_potential_norm
    else
       k_u = 0.0_dp
       k_phi = 0.0_dp
       k_t = 0.0_dp
    end if

    if(.not.all(ieee_is_finite((/h_u,k_u,h_phi,k_phi,h_t,k_t/)))) &
        error stop 'non-finite controlled Love number'
    if(.not.all(ieee_is_finite(residuals))) &
        error stop 'non-finite controlled residual'

    reciprocity_error = surface_gravity*acceleration_norm*h_phi-k_u
    if(max(abs(surface_gravity*acceleration_norm*h_phi),abs(k_u)) > &
       0.0_dp) then
       reciprocity_error = reciprocity_error / &
           max(abs(surface_gravity*acceleration_norm*h_phi),abs(k_u))
    end if

    write(*,'(i0,1x,8(es26.17e3,1x))') &
        degree,h_u,k_u,h_phi,k_phi,h_u+h_phi,k_u+k_phi,h_t,k_t
    write(error_unit,'(a,1x,2(i0,1x),4(es26.17e3,1x))') &
        'diagnostic',degree,matrix%ndim,residuals,reciprocity_error
  end subroutine solve_degree


  subroutine relative_residual(band_matrix,solution,right_hand_side, &
                               bandwidth,residual)
    real(dp), dimension(:,:), intent(in) :: band_matrix
    real(dp), dimension(:), intent(in) :: solution
    real(dp), dimension(:), intent(in) :: right_hand_side
    integer(i4b), intent(in) :: bandwidth
    real(dp), intent(out) :: residual

    real(dp), dimension(:), allocatable :: difference

    allocate(difference(size(right_hand_side)))
    difference = -right_hand_side
    call dsbmv('U',size(solution),bandwidth,1.0_dp,band_matrix, &
               size(band_matrix,1),solution,1,1.0_dp,difference,1)
    residual = sqrt(sum(difference*difference)) / &
               max(1.0_dp,sqrt(sum(right_hand_side*right_hand_side)))
  end subroutine relative_residual

end program controlled_love_numbers
