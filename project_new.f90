module dynamics_module1
        implicit none
        integer::i,j,k,l,Au_nei_list(528,12),step_dynamics,step_equilibrium,step,vib_state,Au_i(528),n_O_Au,n_N_Au,n_Au,traj
        real*8,dimension(396,3)::acc_Au,acc_Au_old,Au_vel
        real*8,dimension(528,3)::Au,Au_eq
        real*8,dimension(3)::N,O,N_vel,O_vel,COM,COM_v,acc_N,acc_N_old,acc_O,acc_O_old
        real*8::r_Au_O,r_Au_N
        integer,dimension(528)::N_Au_list,O_Au_list
        real*8::r_N_O

        real*8::a0,alpha0,b0,beta0,f0,gamma0,r0_N_O,a1,alpha1,b1,beta1,r1_Au_N,d,c,zimage,f1,gamma1,r1_N_O,phi,ea
        real*8::a2,a3,gamma2,gamma3,b2,b3,r_cut
        
        real*8::mass_N,mass_O,red_mass,av,kb,temp,conv,wavenumber_2_j,amu2kg,c_light,wavenumber,mass_Au
        real*8::e_kin,pi,planck_const,total_mass
        
        real*8,dimension(528,528,3,3)::ten_mat
        real*8::l_x,r_x,l_y,r_y,l_z,r_z,box_len_x,box_len_y,box_len_z,force_const,dt1,dt,z_COM
        integer::iflag_write,iflag_writevmd,adsorbed,no_traj
        
        real*8::diab_pot(2,2),gr_vec(1,2),gr_vec_t(2,1),ex_vec(1,2),ex_vec_t(2,1),gr_eig,ex_eig
        real*8::der_pot_N(2,2,3),der_pot_O(2,2,3),der_pot_Au(528,3,2,2),pot_Au_Au(1,1),der_pot_Au_Au(528,3)

contains

subroutine main_sub
        implicit none
        integer::i

        adsorbed=0
        no_traj=100
        !open(7,file='dynamics.xyz')
        call Au_initial_coordinate
        call parameter_value
        call tensor_matrix1
        do traj=1,no_traj
           call gold_equilibrium
           call dynamics
           if(COM(3)<r_cut) adsorbed=adsorbed+1
        enddo
        write(*,*)'for',e_kin/conv,'prob of adsorption',(adsorbed/real(no_traj))
        

end subroutine main_sub


subroutine parameter_value
        implicit none
        pi=4.d0*datan(1.d0)
        amu2kg=1.66053892d-27
        av=6.0221367d23
        conv=1000.d0/av
        kb=8.314d0/av
        mass_O=15.999d0*amu2kg
        mass_N=14.0067d0*amu2kg
        mass_Au=196.96d0*amu2kg
        red_mass=(mass_O*mass_N)/(mass_O+mass_N)
        total_mass=mass_O+mass_N
        planck_const=6.62607d-34
        c_light=2.99792458d10
        wavenumber_2_j=c_light*planck_const

        a0=457095d0*conv
        alpha0=3.7594*1.d10
        b0=30707*conv
        beta0=3.0082*1.d10
        f0=638.5*conv
        gamma0=2.743*1.d10
        r0_N_O=1.15077d-10

        a1=a0
        alpha1=alpha0
        b1=24.056*conv
        beta1=1.9649*1.d10
        r1_Au_N=2.3491*1.d-10
        d=347.22*conv*1.d-10
        c=1.2423*1.d-10
        zimage=1.1536*1.d-10
        f1=495.98*conv
        gamma1=2.4890*1.d10
        r1_N_O=1.2904*1.d-10
        phi=511.37*conv
        ea=-0.67540*conv
        
        a2=11.842*conv
        a3=0.0061803
        gamma2=1.3693*1.d10
        b2=50*conv
        b3=0.0047213
        gamma3=2.0194*1.d10
        r_cut=10.d-10

        l_x=2.9499787804321573E-010
        r_x=3.3924755974969813E-009
        l_y=8.5158552149309472E-011
        r_y=3.2360249816737607E-009
        l_z=-7.230000000000003E-010
        r_z=0.d0
        box_len_x=r_x-l_x
        box_len_y=r_y-l_y
        box_len_z=r_z-l_z


        iflag_write=1
        iflag_writevmd=0
        step_dynamics=010000    !run the dynamics for 10 ps.
        step_equilibrium=05000  !equilibrium dynamics for 5 ps
        temp=300.d0
        dt1=1.d-15              !time step for equilibration
        dt=1.d-15               !time step for running dynamics
        force_const=1550        !unit N/m
        vib_state=0
        e_kin=10.00*conv        !unit kj/mol
        wavenumber=1903.7d0     !unit cm-1
        n_Au=396
        !write(*,*)e_kin/conv
end subroutine parameter_value


subroutine write_output(iflag)
        implicit none
        integer,intent(in)::iflag
        if(iflag==0)then
        endif
        if(iflag==1)then
           write(7,*)'530'
           write(7,*)
           write(7,*)'N',N(:)*1.d10
           write(7,*)'O',O(:)*1.d10
           do i=1,528
              write(7,*)'AU',Au(i,:)*1.d10
           enddo
        endif
end subroutine write_output


subroutine Au_initial_coordinate
        implicit none
        
        Au_eq=0.d0
        open(1,file='528atom.dat')
        do i=1,528
        read(1,*)Au_eq(i,:)
        enddo
        
        Au_nei_list=0
        open(2,file='gold_neigh_pos.dat')
        do i=1,528
        read(2,*)Au_i(i),Au_nei_list(i,:)
        enddo
end subroutine Au_initial_coordinate


subroutine gaussian_random_number(gau_dis)
  implicit none
  real*8 rnd_num,x1,x2,gau_dis
  call random_number(x1)
  call random_number(x2)
  gau_dis=(sqrt(-2.d0*dlog(x1))*dcos(2.d0*pi*x2))

end subroutine gaussian_random_number


subroutine check_pbc(x,y)               !for gold atoms within boundary
        implicit none
        real*8,intent(in)::x(1,3)
        real*8,intent(out)::y(1,3)
        y=x
        if(x(1,1)<l_x) y(1,1)=x(1,1)+box_len_x
        if(x(1,1)>r_x) y(1,1)=x(1,1)-box_len_x
        if(x(1,2)<l_y) y(1,2)=x(1,2)+box_len_y
        if(x(1,2)>r_y) y(1,2)=x(1,2)-box_len_y
        return
end subroutine check_pbc


subroutine check_pbc1(x,y)              !check pbc of r vector
        implicit none
        real*8,intent(in)::x(1,3)
        real*8,intent(out)::y(1,3)
        if(abs(x(1,1))>box_len_x*0.5d0)then
                if(x(1,1)>0.d0) y(1,1)=x(1,1)-box_len_x
                if(x(1,1)<0.d0) y(1,1)=x(1,1)+box_len_X
        else
                y(1,1)=x(1,1)
        endif
        if(abs(x(1,2))>box_len_y*0.5d0)then
                if(x(1,2)>0.d0) y(1,2)=x(1,2)-box_len_y
                if(x(1,2)<0.d0) y(1,2)=x(1,2)+box_len_y
        else
                y(1,2)=x(1,2)
        endif
        y(1,3)=x(1,3)
end subroutine check_pbc1


subroutine tensor_matrix1
        implicit none
        real*8::rij(1,3),t(3,3),t_trans(3,3),t_mat(3,3),ll,a,thetay,thetaz,rij1(1,3)
        real*8::alpha=-4.94d0,beta=17.15d0,gamma4=19.4d0!,ten_mat(528,528,3,3)
        ten_mat(:,:,:,:)=0.d0
        do i=1,528
        res: do j=1,12
                if(Au_nei_list(i,j)==0) cycle
                k=Au_nei_list(i,j)
                rij(1,:)=(Au_eq(i,:)-Au_eq(k,:))
                call check_pbc1(rij,rij1)
                a=dsqrt(sum(rij1*rij1))

                thetaz=dasin(rij1(1,2)/a)
                !thetay=datan(rij1(1,3)/rij1(1,1))
                thetay=dacos(rij1(1,1)/dsqrt(rij1(1,1)**2+rij1(1,3)**2))
                t(1,1)=dcos(thetay)*dcos(thetaz)
                t(1,2)=dcos(thetay)*dsin(thetaz)
                t(1,3)=dsin(thetay)
                t(2,1)=-dsin(thetaz)
                t(2,2)=dcos(thetaz)
                t(2,3)=0.d0
                t(3,1)=-dsin(thetay)*dcos(thetaz)
                t(3,2)=-dsin(thetay)*dsin(thetaz)
                t(3,3)=dcos(thetay)
                t_mat=0.d0
                t_mat(1,1)=beta+gamma4
                t_mat(2,2)=beta-gamma4
                t_mat(3,3)=alpha
                t=transpose(t)
                t_trans=transpose(t)
                ten_mat(i,k,:,:)=matmul(t_trans,matmul(t_mat,t))


             enddo res
         enddo
end subroutine tensor_matrix1


subroutine pbc_distance(x,y,rij)
  implicit none
  real*8,intent(in)::x(3),y(3)
  real*8,intent(out)::rij(3)

  rij=x-y
  if(rij(1)>box_len_x/2.d0) rij(1)=rij(1)-box_len_x
  if(rij(1)<-box_len_x/2.d0) rij(1)=rij(1)+box_len_x

  if(rij(2)>box_len_y/2.d0) rij(2)=rij(2)-box_len_y
  if(rij(2)<-box_len_y/2.d0) rij(2)=rij(2)+box_len_y

end subroutine pbc_distance

subroutine interaction_list
        implicit none
        real*8::dist(3),r
        N_Au_list=0
        n_N_Au=1
        do i=1,528
        call pbc_distance(N(:),Au(i,:),dist(:))
        r=dsqrt(sum(dist(:)*dist(:)))
        if(r<r_cut)then
                N_Au_list(n_N_Au)=i
                n_N_Au=n_N_Au+1
        endif
        enddo
        O_Au_list=0
        n_O_Au=1
        do i=1,528
        call pbc_distance(O(:),Au(i,:),dist(:))
        r=dsqrt(sum(dist(:)*dist(:)))
        if(r<r_cut)then
                O_Au_list(n_O_Au)=i
                n_O_Au=n_O_Au+1
        endif
        enddo

end subroutine interaction_list


subroutine nitrogen_oxygen_int
        implicit none
        real*8::tmp1,dist(3),deriv(3)

        call pbc_distance(N(:),O(:),dist(:))
        r_N_O=dsqrt(sum(dist(:)*dist(:)))
        dist=dist/r_N_O

        tmp1=dexp(-gamma0*(r_N_O-r0_N_O))
        deriv(:)=2.d0*f0*(1.d0-tmp1)*tmp1*gamma0*dist(:)
        diab_pot(1,1)=diab_pot(1,1)+f0*((1.d0-tmp1)**2)
        der_pot_N(1,1,:)=der_pot_N(1,1,:)+deriv(:)
        der_pot_O(1,1,:)=der_pot_O(1,1,:)-deriv(:)

        tmp1=dexp(-gamma1*(r_N_O-r1_N_O))
        deriv(:)=2.d0*f1*(1.d0-tmp1)*gamma1*tmp1*dist(:)
        diab_pot(2,2)=diab_pot(2,2)+f1*((1.d0-tmp1)**2)
        der_pot_N(2,2,:)=der_pot_N(2,2,:)+deriv(:)
        der_pot_O(2,2,:)=der_pot_O(2,2,:)-deriv(:)

end subroutine nitrogen_oxygen_int



subroutine gold_oxygen_int
        implicit none
        real*8::tmp1,dist(3),deriv(3)

        do i=1,n_O_Au
           j=O_Au_list(i)
           call pbc_distance(O(:),Au(j,:),dist(:))
           r_Au_O=dsqrt(sum(dist(:)*dist(:)))
           dist(:)=dist(:)/r_Au_O

           tmp1=dexp(-alpha0*r_Au_O)
           deriv(:)=a0*tmp1*(-alpha0)*dist(:)
           diab_pot(1,1)=diab_pot(1,1)+a0*(tmp1-dexp(-alpha0*r_cut))
           der_pot_O(1,1,:)=der_pot_O(1,1,:)+deriv(:)
           der_pot_Au(j,:,1,1)=der_pot_Au(j,:,1,1)-deriv(:)

           tmp1=dexp(-alpha1*r_Au_O)
           deriv(:)=a1*tmp1*(-alpha1)*dist(:)
           diab_pot(2,2)=diab_pot(2,2)+a1*(tmp1-dexp(-alpha1*r_cut))
           der_pot_O(2,2,:)=der_pot_O(2,2,:)+deriv(:)
           der_pot_Au(j,:,2,2)=der_pot_Au(j,:,2,2)-deriv(:)

           tmp1=dexp(gamma2*r_Au_O)
           deriv(:)=a2*a3*gamma2*tmp1*dist(:)/((1.d0+a3*tmp1)**2)
           diab_pot(1,2)=diab_pot(1,2)-a2*((1.d0/(1.d0+a3*tmp1))-(1.d0/(1.d0+a3*dexp(gamma2*r_cut))))
           der_pot_O(1,2,:)=der_pot_O(1,2,:)+deriv(:)
           der_pot_Au(j,:,1,2)=der_pot_Au(j,:,1,2)-deriv(:)
         enddo
end subroutine gold_oxygen_int



subroutine gold_nitrogen_int
        implicit none
        real*8::tmp1,deriv(3),dist(3),costheta,dcostheta(3)

        do i=1,n_N_Au
           j=N_Au_list(i)
           call pbc_distance(N(:),Au(j,:),dist(:))
           r_Au_N=dsqrt(sum(dist(:)*dist(:)))
           dist=dist/r_Au_N

           tmp1=dexp(-beta0*r_Au_N)
           deriv(:)=b0*tmp1*(-beta0)*dist(:)
           diab_pot(1,1)=diab_pot(1,1)+b0*(tmp1-dexp(-beta0*r_cut))
           der_pot_N(1,1,:)=der_pot_N(1,1,:)+deriv(:)
           der_pot_Au(j,:,1,1)=der_pot_Au(j,:,1,1)-deriv(:)

           tmp1=dexp(-2.d0*beta1*(r_Au_N-r1_Au_N))
           deriv(:)=b1*tmp1*(-2.d0*beta1)*dist(:)
           diab_pot(2,2)=diab_pot(2,2)+b1*(tmp1-dexp(-2.d0*beta1*(r_cut-r1_Au_N)))
           der_pot_N(2,2,:)=der_pot_N(2,2,:)+deriv(:)
           der_pot_Au(j,:,2,2)=der_pot_Au(j,:,2,2)-deriv(:)

           tmp1=dexp(-beta1*(r_Au_N-r1_Au_N))
           deriv(:)=b1*tmp1*(-beta1)*dist(:)
           costheta=(O(3)-N(3))/r_N_O
           dcostheta(:)=-((costheta)*(N(:)-O(:)))/(r_N_O**2)
           dcostheta(3)=dcostheta(3)-1.d0/r_N_O
           diab_pot(2,2)=diab_pot(2,2)-2.d0*b1*(costheta**2)*(tmp1-dexp(-beta1*(r_cut-r1_Au_N)))
           der_pot_N(2,2,:)=der_pot_N(2,2,:)-4.d0*b1*costheta*(tmp1-dexp(-beta1*(r_cut-r1_Au_N)))*dcostheta(:)
           der_pot_O(2,2,:)=der_pot_O(2,2,:)-4.d0*b1*costheta*(tmp1-dexp(-beta1*(r_cut-r1_Au_N)))*(-dcostheta(:))
           der_pot_Au(j,:,2,2)=der_pot_Au(j,:,2,2)+2.d0*(costheta**2)*deriv(:)
           der_pot_N(2,2,:)=der_pot_N(2,2,:)-2.d0*(costheta**2)*deriv(:)

           tmp1=dexp(gamma3*r_Au_N)
           deriv(:)=b2*b3*tmp1*gamma3*dist(:)/((1.d0+b3*tmp1)**2)
           diab_pot(1,2)=diab_pot(1,2)-b2*((1.d0/(1.d0+b3*tmp1))-(1.d0/(1.d0+b3*dexp(gamma3*r_cut))))
           der_pot_N(1,2,:)=der_pot_N(1,2,:)+deriv(:)
           der_pot_Au(j,:,1,2)=der_pot_Au(j,:,1,2)-deriv(:)
           enddo

end subroutine gold_nitrogen_int



subroutine gold_gold_int
        implicit none
        real*8::pot(1,3),dist(1,3),dist_eq(1,3),tensor_mat(3,3),dist_trans(3,1)
        integer::l,m
        der_pot_Au_Au(:,:)=0.d0
        pot_Au_Au=0.d0
        do i=1,528
           res:do j=1,12
               tensor_mat(:,:)=0.d0
               if(Au_nei_list(i,j)==0) cycle
               k=Au_nei_list(i,j)
               call pbc_distance(Au_eq(i,:),Au_eq(k,:),dist_eq(1,:))
               call pbc_distance(Au(i,:),Au(k,:),dist(1,:))
               do l=1,3
               do m=1,3
               tensor_mat(l,m)=ten_mat(i,k,l,m)
               enddo
               enddo
               dist=dist-dist_eq
               dist_trans=transpose(dist)
               pot=matmul(dist,tensor_mat)
               pot_Au_Au=pot_Au_Au+0.5d0*(matmul(pot,dist_trans))
               der_pot_Au_Au(i,:)=der_pot_Au_Au(i,:)+pot(1,:)
               der_pot_Au_Au(k,:)=der_pot_Au_Au(k,:)-pot(1,:)
               enddo res
         enddo

end subroutine gold_gold_int


subroutine image_potential
        implicit none
        real*8::tmp1,tmp2
        
        z_COM=(mass_N*N(3)+mass_O*O(3))/(mass_N+mass_O)
        tmp1=dsqrt(c*c+(z_COM-zimage)**2)
        diab_pot(2,2)=diab_pot(2,2)-d/tmp1
        der_pot_N(2,2,3)=d*(z_COM-zimage)*(mass_N/total_mass)/(tmp1**3)+der_pot_N(2,2,3)
        der_pot_O(2,2,3)=d*(z_COM-zimage)*(mass_O/total_mass)/(tmp1**3)+der_pot_O(2,2,3)
        
end subroutine image_potential


subroutine velocity_verlet
        implicit none
        real*8::surf_temp,multiplier,k_Au
        integer::file_no

        N(:)=N(:)+N_vel(:)*dt+0.5d0*(acc_N(:))*dt**2
        O(:)=O(:)+O_vel(:)*dt+0.5d0*(acc_O(:))*dt**2
        do i=1,396
        Au(i,:)=Au(i,:)+Au_vel(i,:)*dt+0.5d0*(acc_Au(i,:))*dt**2
        enddo
        
        file_no=16+traj
        COM(:)=(mass_N*N(:)+mass_O*O(:))/total_mass
        !write(file_no,*)step*dt*1.d12,COM(3)*1.d10,N(3)*1.d10,O(3)*1.d10
         
        call check_pbc(N(:),N(:))
        call check_pbc(O(:),O(:))
        do i=1,396
        call check_pbc(Au(i,:),Au(i,:))
        enddo
        !call write_output(1)

        acc_N_old=acc_N
        acc_O_old=acc_O
        acc_Au_old=acc_Au

        call calculate_acc

        N_vel=N_vel+0.5d0*(acc_N+acc_N_old)*dt
        O_vel=O_vel+0.5d0*(acc_O+acc_O_old)*dt
        
        k_Au=0.d0
        do i=1,396
                Au_vel(i,:)=Au_vel(i,:)+0.5d0*(acc_Au(i,:)+acc_Au_old(i,:))*dt
        enddo
        
        call gold_surf_temp(multiplier)
        Au_vel=Au_vel*multiplier
end subroutine velocity_verlet

subroutine total_energy
        implicit none
        real*8::energy,k_Au,ke_NO,vel(3),ke_vib,ke_rot,multiplier
        
        energy=0.d0
        energy=pot_Au_Au(1,1)
        k_Au=0.d0
        call gold_surf_temp(multiplier)
        ke_NO=0.5d0*mass_N*sum(N_vel(:)*N_vel(:))+0.5d0*mass_O*sum(O_vel(:)*O_vel(:))
        energy=energy+k_Au+ke_NO

        r_N_O=dsqrt(sum((N(:)-O(:))**2))
        vel=N_vel-O_vel
        ke_vib=0.5*red_mass*(sum(vel*vel))**2
        ke_rot=0.5d0*red_mass*r_N_O**2*((sum(vel*vel)/r_N_O)**2)
        
end subroutine


subroutine gold_surf_temp(multiplier)
        implicit none
        real*8::k_Au,surf_temp
        real*8,intent(out)::multiplier
        k_Au=0.d0
        do i=1,396
        do j=1,3
        k_Au=k_Au+0.5d0*mass_Au*((Au_vel(i,j))**2)
        enddo
        enddo

        surf_temp=k_Au*2.d0/(3.d0*kb*396)
        multiplier=dsqrt(temp/surf_temp)
        !write(12,*)step,surf_temp,k_Au/conv,pot_Au_Au/conv,multiplier
end subroutine 


subroutine dynamics
        implicit none
        real*8::multiplier
        
        step=0
        call set_NO2
        call calculate_acc

        do step=1,step_dynamics
                call velocity_verlet
        enddo
        
        COM(3)=(mass_N*N(3)+mass_O*O(3))/(mass_N+mass_O)
end subroutine dynamics


subroutine gold_equilibrium
        implicit none
        real*8::k_Au,multiplier

        Au(:,:)=Au_eq(:,:)
        call gold_eq_acc
        call gold_init_vel
        call gold_surf_temp(multiplier)
        Au_vel=Au_vel*multiplier
        acc_Au_old(:,:)=0.d0
        do step=1,step_equilibrium
             
             do i=1,396
             do j=1,3
             Au(i,j)=Au(i,j)+Au_vel(i,j)*dt1+0.5d0*(acc_Au(i,j))*dt1**2
             enddo
             enddo
             do i=1,396
             call check_pbc(Au(i,:),Au(i,:))
             enddo
             
             acc_Au_old(:,:)=acc_Au(:,:)
             call gold_eq_acc
             k_Au=0.d0
             do i=1,396
             do j=1,3
             Au_vel(i,j)=Au_vel(i,j)+0.5d0*(acc_Au_old(i,j)+acc_Au(i,j))*dt1
             !k_Au=k_Au+0.5d0*mass_Au*(Au_vel(i,j)**2)
             enddo
             enddo
             !write(12,*)step,k_Au/conv,pot_Au_Au/conv
         enddo

end subroutine gold_equilibrium


subroutine gold_vel_scalling
        implicit none
        real*8::total_vel(3)

        total_vel=0.d0
        do i=1,396
        total_vel(:)=total_vel(:)+Au_vel(i,:)
        enddo
        total_vel(:)=total_vel(:)/396
        do i=1,396
        Au_vel(i,:)=Au_vel(i,:)-total_vel(:)
        enddo

end subroutine gold_vel_scalling

subroutine gold_init_vel
        implicit none
        real*8::sig_p,rnd
        
        Au_vel=0.d0
        sig_p=dsqrt(kb*temp/mass_Au)
        do i=1,396
           do j=1,3
           call gaussian_random_number(rnd)
           Au_vel(i,j)=(rnd*sig_p)!+sig_p
           enddo
        enddo

        call gold_vel_scalling
end subroutine gold_init_vel

subroutine gold_eq_acc
        implicit none

        call gold_gold_int
        do i=1,396
        acc_Au(i,:)=-(der_pot_Au_Au(i,:))/mass_Au
        enddo

end subroutine gold_eq_acc


subroutine mat_diag(h,gr_eig,ex_eig,gr_vec,ex_vec)
        implicit none
        real*8,dimension(2,2),intent(in)::h
        real*8,intent(out)::gr_eig,ex_eig
        real*8,dimension(1,2),intent(out)::gr_vec,ex_vec
        real*8::tmp
        
        gr_eig=0.5d0*((h(1,1)+h(2,2))-dsqrt((h(1,1)-h(2,2))**2+4.d0*(h(1,2))**2))
        ex_eig=0.5d0*((h(1,1)+h(2,2))+dsqrt((h(1,1)-h(2,2))**2+4.d0*(h(1,2))**2))

        tmp=(h(1,2)**2)+(h(1,1)-gr_eig)**2
        gr_vec(1,1)=h(1,2)/dsqrt(tmp)
        gr_vec(1,2)=-(h(1,1)-gr_eig)/dsqrt(tmp)
        !write(13,*)step,h(1,1)/conv,h(1,2)/conv,gr_eig/conv

        tmp=(h(1,2)**2)+(h(1,1)-ex_eig)**2
        ex_vec(1,1)=h(1,2)/dsqrt(tmp)
        ex_vec(1,2)=(ex_eig-h(1,1))/dsqrt(tmp)

end subroutine mat_diag


subroutine calculate_acc
        implicit none
        real*8::c11,c12
        
        diab_pot=0.d0
        der_pot_N=0.d0
        der_pot_O=0.d0
        der_pot_Au=0.d0
        acc_N=0.d0
        acc_O=0.d0
        acc_Au=0.d0
        
        call interaction_list
        call nitrogen_oxygen_int
        call gold_nitrogen_int
        call gold_oxygen_int
        call image_potential
        call gold_gold_int

        diab_pot(2,1)=diab_pot(1,2)
        diab_pot(2,2)=diab_pot(2,2)+phi-ea
        call mat_diag(diab_pot,gr_eig,ex_eig,gr_vec,ex_vec)
        c11=gr_vec(1,1)
        c12=gr_vec(1,2)
        !write(12,*)step,c11,c12,gr_eig/conv
 
        der_pot_N(2,1,:)=der_pot_N(1,2,:)
        der_pot_O(2,1,:)=der_pot_O(1,2,:)
        der_pot_Au(:,:,2,1)=der_pot_Au(:,:,1,2)


        acc_N(:)=-((der_pot_N(1,1,:)*c11**2)+(2.d0*der_pot_N(1,2,:)*(c11*c12))+(der_pot_N(2,2,:)*c12**2))/mass_N
        acc_O(:)=-((der_pot_O(1,1,:)*c11**2)+(2.d0*der_pot_O(1,2,:)*(c11*c12))+(der_pot_O(2,2,:)*c12**2))/mass_O
        do i=1,396
           acc_Au(i,:)=-((der_pot_Au(i,:,1,1)*(c11**2))+(2.d0*der_pot_Au(i,:,1,2)*(c11*c12)))/mass_Au+acc_Au(i,:)
           acc_Au(i,:)=-((der_pot_Au(i,:,2,2)*(c12**2))+der_pot_Au_Au(i,:))/mass_Au+acc_Au(i,:)
        enddo
end subroutine calculate_acc
        

subroutine set_NO1              !NO is perpendicular to the surface. vibration only in z direction
        implicit none
        real*8::rel_x,rel_v,theta,phi,rnd,omega0,vib_energy,r_N_z,r_O_z,rand_ori
        
        call random_number(rnd)
        COM(1)=l_x+box_len_x*rnd*0.33d0+box_len_x*0.33d0
        call random_number(rnd)
        COM(2)=l_y+box_len_y*rnd*0.33d0+box_len_y*0.33d0
        COM(3)=10d-10
        
        omega0=dsqrt(2.d0*f0*gamma0**2/red_mass)
        !vib_energy=(real(vib_state)+0.5d0)*wavenumber*wavenumber_2_j
        vib_energy=0.5*planck_const*omega0/(2.d0*pi)
!        write(*,*)'vib_energy',vib_energy/conv
        call random_number(rnd)
        theta=pi*rnd*2.d0
        rel_x=dsqrt(2.d0*vib_energy/(red_mass*omega0**2))*dcos(theta)
        rel_v=dsqrt(2.d0*vib_energy/red_mass)*dsin(theta)
        write(9,*)'rel_x',rel_x,'rel_v',rel_v
        call random_number(rnd)
        if(rnd<0.5d0) rand_ori=1.d0
        if(rnd>0.5d0) rand_ori=-1.d0
        r_N_z=r0_N_O*mass_O/total_mass
        r_O_z=r0_N_O*mass_N/total_mass

        N(1)=COM(1)
        N(2)=COM(2)
        N(3)=COM(3)-(r0_N_O+rel_x)*mass_O/total_mass

        O(1)=COM(1)
        O(2)=COM(2)
        O(3)=COM(3)+(r0_N_O+rel_x)*mass_N/total_mass

        N_vel(1)=0.d0
        N_vel(2)=0.d0
        N_vel(3)=-rel_v*0.5d0
        O_vel(1)=0.d0
        O_vel(2)=0.d0
        O_vel(3)=+rel_v*0.5d0

        N_vel(3)=N_vel(3)-dsqrt(2.d0*e_kin/total_mass)
        O_vel(3)=O_vel(3)-dsqrt(2.d0*e_kin/total_mass)
end subroutine set_NO1

subroutine set_NO2                      !any random NO orientation. random theta and random phi.
        implicit none
        real*8::rel_x,rel_v,theta,phi,rnd,omega0,vib_energy

        call random_number(rnd)
        COM(1)=l_x+box_len_x*rnd!*0.33d0+box_len_x*0.33d0
        call random_number(rnd)
        COM(2)=l_y+box_len_y*rnd!*0.33d0+box_len_y*0.33d0
        COM(3)=11.d-10

        omega0=dsqrt(2.d0*f0*gamma0**2/red_mass)
        !omega0=dsqrt(force_const/red_mass)
        !vib_energy=(real(vib_state)+0.5d0)*wavenumber*wavenumber_2_j
        vib_energy=(0.5+real(vib_state))*planck_const*omega0/(2.d0*pi)

        call random_number(rnd)
        theta=2.d0*pi*rnd
        rel_x=dsqrt(2.d0*vib_energy/(red_mass*omega0**2))*dcos(theta)
        rel_v=-dsqrt(2.d0*vib_energy/red_mass)*dsin(theta)

        !call random_sphere(theta,phi)          this does not provide complete randomness
        call random_number(rnd)
        theta=pi*rnd
        call random_number(rnd)
        phi=2.d0*pi*rnd
        
        N(1)=COM(1)-(r0_N_O+rel_x)*dsin(theta)*dcos(phi)*mass_O/(mass_O+mass_N)
        N(2)=COM(2)-(r0_N_O+rel_x)*dsin(theta)*dsin(phi)*mass_O/(mass_O+mass_N)
        N(3)=COM(3)-(r0_N_O+rel_x)*dcos(theta)*mass_O/(mass_O+mass_N)

        O(1)=COM(1)+(r0_N_O+rel_x)*dsin(theta)*dcos(phi)*mass_N/(mass_O+mass_N)
        O(2)=COM(2)+(r0_N_O+rel_x)*dsin(theta)*dsin(phi)*mass_N/(mass_O+mass_N)
        O(3)=COM(3)+(r0_N_O+rel_x)*dcos(theta)*mass_N/(mass_O+mass_N)

        N_vel(1)=rel_v*dcos(phi)*dsin(theta)*0.5d0*mass_O/total_mass
        N_vel(2)=rel_v*dsin(phi)*dsin(theta)*0.5d0*mass_O/total_mass
        N_vel(3)=-dsqrt(2*e_kin/total_mass)+rel_v*dcos(theta)*0.5d0*mass_O/total_mass

        O_vel(1)=-rel_v*dcos(phi)*dsin(theta)*0.5d0*mass_N/total_mass
        O_vel(2)=-rel_v*dsin(phi)*dsin(theta)*0.5d0*mass_N/total_mass
        O_vel(3)=-dsqrt(2*e_kin/total_mass)-rel_v*dcos(theta)*0.5d0*mass_N/total_mass


end subroutine

subroutine random_sphere(theta,phi)
  implicit none
  real*8,intent(out)::theta,phi
  real*8 vec(3),rnd1,rnd2,rndsq,r
  real*8 pi

  pi=dacos(-1.d0)

  rndsq=2.d0
  do while(rndsq>1.d0)
    call random_number(rnd1)
    call random_number(rnd2)
    rnd1=1.d0-2*rnd1
    rnd2=1.d0-2*rnd2
    rndsq=rnd1**2+rnd2**2
  enddo
  vec(1)=2*rnd1*sqrt(1-rndsq)
  vec(2)=2*rnd2*sqrt(1-rndsq)
  vec(3)=1-2*rndsq

  r=sqrt(sum(vec*vec))
  theta=dacos(vec(3)/r)
  phi=datan(vec(2)/vec(1))
  if(vec(1)>0.d0.and.vec(2)<0.d0)phi=phi+2*pi
  if(vec(1)<0.d0)phi=phi+pi

end subroutine random_sphere

end module dynamics_module1

