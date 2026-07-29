program export_isotropic_prem

  use, intrinsic :: iso_fortran_env, only : error_unit
  use module_constants
  use module_physical_constants, only : density_norm, length_norm, &
                                        velocity_norm
  use module_PREM
  use module_spherical_model
  implicit none

  character(len=512) :: argument,dspec_file,gia3d_file
  integer(i4b) :: dspec_unit,gia3d_unit,ilayer,interval,isection
  integer(i4b) :: intervals,layer_count,knot_count
  real(dp) :: density,kappa,maximum_knot_spacing,mu
  real(dp) :: radius,radius_fraction,p_wave_speed,s_wave_speed
  type(spherical_model) :: model

  if(command_argument_count() /= 3) then
     error stop 'usage: export_isotropic_prem spacing_m dspec_file gia3d_file'
  end if
  call get_command_argument(1,argument)
  read(argument,*) maximum_knot_spacing
  call get_command_argument(2,dspec_file)
  call get_command_argument(3,gia3d_file)
  if(maximum_knot_spacing <= 0.0_dp) &
      error stop 'maximum knot spacing must be positive'

  model = elastic_PREM(.false.)
  layer_count = 0
  knot_count = 0
  do isection = 1,model%nsections
     do ilayer = 1,model%section(isection)%nlayers
        associate(layer => model%section(isection)%layer(ilayer))
          intervals = max(3_i4b,ceiling( &
              (layer%r2-layer%r1)*length_norm/maximum_knot_spacing))
          knot_count = knot_count+intervals+1
          layer_count = layer_count+1
        end associate
     end do
  end do
  if(layer_count /= 12) error stop 'ocean-free PREM must have 12 layers'

  open(newunit=dspec_unit,file=trim(dspec_file),status='replace')
  open(newunit=gia3d_unit,file=trim(gia3d_file),status='replace')
  write(dspec_unit,'(a)') &
      'isotropized analytic no-ocean PREM exported from pinned gia3D'
  write(dspec_unit,'(a)') '1 -1.0 1 2'
  write(dspec_unit,'(i0,a)') knot_count,' 0 0'

  do isection = 1,model%nsections
     do ilayer = 1,model%section(isection)%nlayers
        associate(layer => model%section(isection)%layer(ilayer))
          intervals = max(3_i4b,ceiling( &
              (layer%r2-layer%r1)*length_norm/maximum_knot_spacing))
          do interval = 0,intervals
             radius_fraction = real(interval,dp)/real(intervals,dp)
             radius = layer%r1+radius_fraction*(layer%r2-layer%r1)
             density = layer%rho(radius)

             select type(layer)
             class is(spherical_solid_elastic_layer)
                kappa = (layer%C(radius)+4.0_dp*layer%A(radius) - &
                    4.0_dp*layer%N(radius)+4.0_dp*layer%F(radius))/9.0_dp
                mu = (layer%C(radius)+layer%A(radius) + &
                    6.0_dp*layer%L(radius)+5.0_dp*layer%N(radius) - &
                    2.0_dp*layer%F(radius))/15.0_dp
                p_wave_speed = sqrt((kappa+4.0_dp*mu/3.0_dp)/density)
                s_wave_speed = sqrt(mu/density)
             class is(spherical_fluid_elastic_layer)
                p_wave_speed = sqrt(layer%kappa(radius)/density)
                s_wave_speed = 0.0_dp
             class default
                error stop 'unsupported PREM layer in isotropic exporter'
             end select

             radius = radius*length_norm
             density = density*density_norm
             p_wave_speed = p_wave_speed*velocity_norm
             s_wave_speed = s_wave_speed*velocity_norm
             write(gia3d_unit,'(4(es26.17e3,1x))') &
                 radius,density,p_wave_speed,s_wave_speed
             write(dspec_unit,'(9(es26.17e3,1x))') &
                 radius,density,p_wave_speed,s_wave_speed, &
                 1000.0_dp,600.0_dp,p_wave_speed,s_wave_speed,1.0_dp
          end do
        end associate
     end do
  end do
  close(dspec_unit)
  close(gia3d_unit)

  write(error_unit,'(a,1x,2(i0,1x),2(es26.17e3,1x))') &
      'exported',layer_count,knot_count,maximum_knot_spacing, &
      model%r2*length_norm

end program export_isotropic_prem
