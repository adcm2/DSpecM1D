program direct_love_numbers

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : error_unit
  use module_constants
  use module_physical_constants, only : gravitational_potential_norm, &
                                        length_norm, load_norm
  use module_PREM
  use module_spherical_model
  use module_spherical_mesh
  use module_matrix
  use module_force
  implicit none

  character(len=32) :: argument
  integer(i4b) :: argument_index, maximum_degree, number_of_degrees
  integer(i4b), dimension(:), allocatable :: degrees
  real(dp) :: maximum_radial_step
  type(spherical_model) :: model
  type(spherical_model_mesh) :: mesh

  number_of_degrees = command_argument_count()
  if(number_of_degrees == 0) error stop 'provide at least one degree'
  allocate(degrees(number_of_degrees))

  do argument_index = 1,number_of_degrees
     call get_command_argument(argument_index,argument)
     read(argument,*) degrees(argument_index)
     if(degrees(argument_index) < 1) error stop 'degrees must be positive'
  end do

  maximum_degree = maxval(degrees)
  model = elastic_PREM(.false.)
  maximum_radial_step = 0.1_dp*model%r2/real(maximum_degree+1,dp)
  mesh = spherical_mesh(5,model,maximum_radial_step)

  do argument_index = 1,number_of_degrees
     call solve_degree(mesh,degrees(argument_index))
  end do

contains

  subroutine solve_degree(mesh,degree)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: degree

    integer(i4b) :: info, surface_element, surface_layer, surface_node
    integer(i4b) :: surface_potential, surface_radial
    real(dp) :: h_load, h_phi, h_t, h_u
    real(dp) :: k_load, k_phi, k_t, k_u
    real(dp) :: stock_h, stock_h_difference, stock_h_relative
    real(dp) :: stock_k, stock_k_difference, stock_k_relative
    real(dp), dimension(3) :: residuals
    real(dp) :: stock_residual, surface_gravity, surface_radius
    real(dp), dimension(:,:), allocatable :: original_matrix
    real(dp), dimension(:,:), allocatable :: right_hand_sides
    real(dp), dimension(:,:), allocatable :: right_hand_side_copy
    real(dp), dimension(:,:), allocatable :: stock_right_hand_side
    real(dp), dimension(:,:), allocatable :: stock_right_hand_side_copy
    type(radial_matrix) :: matrix

    matrix = build_spheroidal_matrix(mesh,degree,.false.)
    allocate(original_matrix(matrix%ldab,matrix%ndim))
    original_matrix = matrix%a

    allocate(right_hand_sides(matrix%ndim,3))
    allocate(right_hand_side_copy(matrix%ndim,3))
    allocate(stock_right_hand_side(matrix%ndim,1))
    allocate(stock_right_hand_side_copy(matrix%ndim,1))
    right_hand_sides = 0.0_dp
    stock_right_hand_side = 0.0_dp

    surface_layer = mesh%section(mesh%nsections)%nlayers
    associate(layer => mesh%section(mesh%nsections)%layer(surface_layer), &
              ibool => matrix%ibool%section(mesh%nsections)%layer(surface_layer))
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

        call force_for_unit_harmonic_load( &
            layer,ibool,degree,stock_right_hand_side)
      class default
        error stop 'gia3D reference requires a solid surface'
      end select
    end associate

    call force_for_unit_harmonic_tide( &
        mesh,matrix%ibool,degree,right_hand_sides(:,3:3))
    right_hand_side_copy = right_hand_sides
    stock_right_hand_side_copy = stock_right_hand_side

    call dpbtrf('U',matrix%ndim,matrix%kd,matrix%a,matrix%ldab,info)
    if(info /= 0) error stop 'gia3D matrix factorisation failed'

    call dpbtrs('U',matrix%ndim,matrix%kd,3,matrix%a,matrix%ldab, &
                right_hand_sides,matrix%ndim,info)
    if(info /= 0) error stop 'gia3D three-column solve failed'
    call dpbtrs('U',matrix%ndim,matrix%kd,1,matrix%a,matrix%ldab, &
                stock_right_hand_side,matrix%ndim,info)
    if(info /= 0) error stop 'gia3D stock loading solve failed'

    call relative_residual(original_matrix,right_hand_sides(:,1), &
                           right_hand_side_copy(:,1),matrix%kd,residuals(1))
    call relative_residual(original_matrix,right_hand_sides(:,2), &
                           right_hand_side_copy(:,2),matrix%kd,residuals(2))
    call relative_residual(original_matrix,right_hand_sides(:,3), &
                           right_hand_side_copy(:,3),matrix%kd,residuals(3))
    call relative_residual(original_matrix,stock_right_hand_side(:,1), &
                           stock_right_hand_side_copy(:,1),matrix%kd, &
                           stock_residual)

    h_u = right_hand_sides(surface_radial,1)*length_norm/load_norm
    h_phi = right_hand_sides(surface_radial,2)*length_norm/load_norm
    h_t = right_hand_sides(surface_radial,3)*length_norm
    stock_h = stock_right_hand_side(surface_radial,1)*length_norm/load_norm

    if(degree > 1) then
       k_u = right_hand_sides(surface_potential,1) * &
             gravitational_potential_norm/load_norm
       k_phi = right_hand_sides(surface_potential,2) * &
               gravitational_potential_norm/load_norm
       k_t = right_hand_sides(surface_potential,3) * &
             gravitational_potential_norm
       stock_k = stock_right_hand_side(surface_potential,1) * &
                 gravitational_potential_norm/load_norm
    else
       k_u = 0.0_dp
       k_phi = 0.0_dp
       k_t = 0.0_dp
       stock_k = 0.0_dp
    end if

    h_load = h_u+h_phi
    k_load = k_u+k_phi
    stock_h_difference = h_load-stock_h
    stock_k_difference = k_load-stock_k
    stock_h_relative = relative_error(h_load,stock_h)
    stock_k_relative = relative_error(k_load,stock_k)

    if(.not.all(ieee_is_finite((/h_u,k_u,h_phi,k_phi,h_load,k_load,h_t,k_t/)))) &
        error stop 'non-finite Love number'
    if(.not.all(ieee_is_finite((/stock_h,stock_k,stock_h_difference, &
                                 stock_k_difference,stock_h_relative, &
                                 stock_k_relative/)))) &
        error stop 'non-finite stock loading comparison'
    if(.not.all(ieee_is_finite(residuals)) .or. &
       .not.ieee_is_finite(stock_residual)) error stop 'non-finite residual'

    write(*,'(i0,1x,8(es26.17e3,1x))') &
        degree,h_u,k_u,h_phi,k_phi,h_load,k_load,h_t,k_t
    write(error_unit,'(a,1x,i0,1x,8(es26.17e3,1x))') &
        'diagnostic',degree,stock_h_difference,stock_h_relative, &
        stock_k_difference,stock_k_relative,residuals,stock_residual
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


  real(dp) function relative_error(value,reference)
    real(dp), intent(in) :: value,reference

    if(reference /= 0.0_dp) then
       relative_error = abs((value-reference)/reference)
    else if(value == reference) then
       relative_error = 0.0_dp
    else
       relative_error = huge(1.0_dp)
    end if
  end function relative_error

end program direct_love_numbers
