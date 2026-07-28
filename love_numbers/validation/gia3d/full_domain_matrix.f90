module full_domain_matrix

  use module_constants
  use module_error
  use module_matrix
  use module_spherical_mesh
  implicit none

contains

  function build_spheroidal_matrix_from_radius( &
      mesh,degree,starting_radius,factor) result(matrix)
    type(spherical_model_mesh), intent(in) :: mesh
    integer(i4b), intent(in) :: degree
    real(dp), intent(in) :: starting_radius
    logical, intent(in), optional :: factor
    type(radial_matrix) :: matrix

    integer(i4b) :: info,isection,ilayer
    integer(i4b) :: first_section,first_layer
    logical :: factor_matrix

    call check(degree > 0,'build_spheroidal_matrix_from_radius', &
               'degree must be positive')

    factor_matrix = .true.
    if(present(factor)) factor_matrix = factor

    matrix%l = degree
    matrix%ibool = build_boolean_spheroidal(mesh,starting_radius)
    if(degree == 1) matrix%ibool%ndim = matrix%ibool%ndim-1
    matrix%ndim = matrix%ibool%ndim
    matrix%kd = 3*(matrix%ibool%ngll-1)+2
    matrix%ldab = matrix%kd+1
    allocate(matrix%a(matrix%ldab,matrix%ndim))
    matrix%a = 0.0_dp

    first_section = matrix%ibool%isection1
    do isection = first_section,mesh%nsections
       first_layer = matrix%ibool%section(isection)%ilayer1
       do ilayer = first_layer,mesh%section(isection)%nlayers
          associate(layer => mesh%section(isection)%layer(ilayer), &
                    ibool => matrix%ibool%section(isection)%layer(ilayer))
            select type(layer)
            class is(spherical_solid_elastic_layer_mesh)
               call build_spheroidal_matrix_solid_layer( &
                   layer,ibool,degree,matrix%a)
            class is(spherical_fluid_elastic_layer_mesh)
               call build_spheroidal_matrix_fluid_layer( &
                   layer,ibool,degree,matrix%a)
            class default
               error stop 'unsupported layer in full-domain matrix'
            end select
          end associate
       end do
    end do

    if(factor_matrix) then
       call dpbtrf('U',matrix%ndim,matrix%kd,matrix%a, &
                   matrix%ldab,info)
       call check(info == 0,'build_spheroidal_matrix_from_radius', &
                  'matrix factorisation failed')
       matrix%factorised = .true.
    end if
  end function build_spheroidal_matrix_from_radius

end module full_domain_matrix
