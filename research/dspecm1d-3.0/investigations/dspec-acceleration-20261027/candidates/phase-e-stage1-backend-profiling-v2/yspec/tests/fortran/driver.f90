program driver
  ! Prints the responses of the radial solvers to unit jumps, for a list of
  ! degrees and frequencies read from standard input:
  !   model_file
  !   ats igrav fic_unused
  !   source_depth_km receiver_depth_km tout_min
  !   n
  !   l f_mHz     (n lines)
  use module_kinds
  use module_model
  use module_solver
  implicit none
  character(len=256) :: model_file
  integer(i4b) :: ats,igrav,n,i,l,is,ir,j,dummy
  real(dp) :: rs,rr,tout,f,ep
  complex(dpc) :: w
  complex(dpc), dimension(2,2) :: s2,ytor,yrad,prad
  complex(dpc), dimension(2) :: sr0
  real(dp), dimension(6) :: mm = 0
  complex(dpc), dimension(4,4) :: s4,y4
  complex(dpc), dimension(6,4) :: s6,y6

  call set_parameters
  read(*,'(a)') model_file
  read(*,*) ats,igrav,dummy
  read(*,*) rs,rr,tout
  open(7,file=trim(model_file),action='read')
  call read_model(7)
  close(7)
  call set_ats_switch(ats)
  rs = r(nknot)-rs*1000.0_dp/r_norm
  rr = r(nknot)-rr*1000.0_dp/r_norm
  is = find_radius_index(rs,-1)
  ir = find_radius_index(rr,-1)
  ep = 5.0_dp/(tout*60.0_dp/t_norm)
  s2 = 0; s2(1,1) = 1; s2(2,2) = 1
  s4 = 0; do j=1,4; s4(j,j) = 1; end do
  s6 = 0; s6(1,1) = 1; s6(2,2) = 1; s6(4,3) = 1; s6(5,4) = 1
  write(*,'(a,es24.16)') '# t_norm ',t_norm
  write(*,'(a,es24.16)') '# ep ',ep
  read(*,*) n
  do i = 1,n
     read(*,*) l,f
     w = twopi_d*f/1000.0_dp*t_norm-ii*ep
     if(l == 0) then
        ! sets the attenuation switch before the solver uses it
        call source_vector_rad(w,is,rs,mm,sr0)
        if(igrav == 0) then
           call rad_solver_ng(w,rr,rs,ir,is,s2,yrad)
        else
           call rad_solver_wg(w,rr,rs,ir,is,s2,yrad,prad)
        end if
        write(*,'(a,i6,es24.16,4es24.16)') 'R',l,f,yrad(1,1),yrad(1,2)
     else
        call tor_solver(l,w,rr,rs,ir,is,s2,ytor)
        write(*,'(a,i6,es24.16,4es24.16)') 'T',l,f,ytor(1,1),ytor(1,2)
        if(fic) then
           if(igrav == 0) call sph_solver_ng_fic(l,w,rr,rs,ir,is,s4,y4)
           if(igrav == 1) call sph_solver_cow_fic(l,w,rr,rs,ir,is,s4,y4)
           if(igrav == 2) call sph_solver_fic(l,w,rr,rs,ir,is,s6,y6)
        else
           if(igrav == 0) call sph_solver_ng(l,w,rr,rs,ir,is,s4,y4)
           if(igrav == 1) call sph_solver_cow(l,w,rr,rs,ir,is,s4,y4)
           if(igrav == 2) call sph_solver(l,w,rr,rs,ir,is,s6,y6)
        end if
        if(igrav /= 2) then
           y6 = 0
           y6(1,:) = y4(1,:); y6(2,:) = y4(2,:)
        end if
        do j = 1,4
           write(*,'(a,i6,es24.16,i3,6es24.16)') 'S',l,f,j,y6(1,j),y6(2,j),y6(3,j)
        end do
     end if
  end do
end program driver
