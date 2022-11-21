module sim_para
    implicit none
    integer           :: zzz=-4280145,zzzz=77777
    integer,parameter :: lx=20,ly=20,lz=20,n_part=2197
    integer,parameter :: kb_T=1,niter=20000,mass=1
    real*8            :: pos(3*n_part),vel(3*n_part),force(3*n_part),acc(3*n_part),new_acc(3*n_part)

    real*8,parameter  :: rc=2.50d0,rs=2.0d0*rc,sigma=1.0d0
    real*8,parameter  :: sigma6=sigma**6,eps=4.0d0,sigma12=sigma**12
    real*8,parameter  :: fc=eps*((12.0d0*sigma12/(rc**13)) -(6.0d0*sigma6/(rc**7) ))
    real*8,parameter  :: ufc=fc*rc+eps*(((sigma/rc)**12) - ((sigma/rc)**6))

    real*8            :: avr_vel_x,avr_vel_y,avr_vel_z
    integer           :: i,j,k,c
    real*8            :: x1,y1,z1,x2,y2,z2,x,y,z,dx,dy,dz,r,lj,pot_energy,ke,gy,gz,gx,delta_t=0.005
    real*8            :: lj_force
    real*8            :: new_pot_energy,new_force(3*n_part)
    real*8            :: invr,ir2,ir6
    real*8            :: llx,lly,llz,theoryke,scalef
endmodule sim_para

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MAIN PROGRAM
program md 
    use sim_para
    implicit none
    real*8 :: ran1,zz1

    zz1=ran1(zzz)
    llx=dfloat(lx);  lly=dfloat(ly);  llz=dfloat(lz)

    call posinit  !INITIALIZE POSITION
    call vel_init !INTIALIZE VEL AND FORCE
    c=0
!         open(90,file='dist.dat',status='unknown')
          open(82,file='energy.dat',status='unknown',form='formatted')
          open(45,file='pot_energy.dat',status='unknown',form='formatted')
          open(47,file='ke_energy.dat',status='unknown',form='formatted')
          
    
    do k=1,niter
        if(mod(k,100)==0) write(*,*) k 
        call pos_update
        write(82,'(4g25.15)') dfloat(k), (new_pot_energy)/dfloat(n_part), ke/dfloat(n_part), (ke+new_pot_energy)/dfloat(n_part)
        write(45,*) dfloat(k), (new_pot_energy)/dfloat(n_part)
        write(47,*) dfloat(k), ke/dfloat(n_part)
    end do 

end program md 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!POSITION & VELOCITY UPDATE
subroutine pos_update
    use sim_para
    implicit none
    real*8 :: zz1,ran1
    real*8 :: dt2by2

! UPDATING POSITION

    dt2by2 = 0.50d0*delta_t**2

    do i=1,n_part
        pos(3*i-2)= pos(3*i-2) +vel(3*i-2)*delta_t + dt2by2*force(3*i-2)
        pos(3*i-1)= pos(3*i-1) +vel(3*i-1)*delta_t + dt2by2*force(3*i-1)
        pos(3*i)= pos(3*i) +vel(3*i)*delta_t + dt2by2*force(3*i)

        pos(3*i-2) = modulo(pos(3*i-2),llx)  !PBC
        pos(3*i-1) = modulo(pos(3*i-1),lly)
        pos(3*i) = modulo(pos(3*i),llz)
    end do 


! CALCULATION OF FORCES

    new_force=0.0d0; new_pot_energy=0.0d0
    acc=force/dfloat(mass)

    do i=1,n_part-1
        x1=pos(3*i-2); y1=pos(3*i-1); z1=pos(3*i)

        do j=i+1,n_part
            c=c+1
            x2=pos(3*j-2); y2=pos(3*j-1); z2=pos(3*j)

            x=x1-x2; y=y1-y2; z=z1-z2 

            if(abs(x) .ge. (dfloat(lx)/2.0d0)) x=(dfloat(lx) - abs(x))*((-1.0d0*x)/abs(x))
            if(abs(y) .ge. (dfloat(ly)/2.0d0)) y=(dfloat(ly) - abs(y))*((-1.0d0*y)/abs(y))
            if(abs(z) .ge. (dfloat(lz)/2.0d0)) z=(dfloat(lz) - abs(z))*((-1.0d0*z)/abs(z))

            r=dsqrt(x*x+y*y+z*z)
            if (r <=rc) then 

                lj= eps*((sigma/r)**12-(sigma/r)**6)-ufc+fc*r
                new_pot_energy = new_pot_energy+lj 
                lj_force= eps*((12.0d0*((sigma12)/(r)**13)) - (6.0d0*((sigma6)/(r)**7))) - fc 

               ! write(90,*) new_force(3*i-2),new_force(3*i-1),new_force(3*i) !Generates around 50gb of data files
                new_force(3*i-2) = new_force(3*i-2) + (lj_force)*(x/r)
                new_force(3*i-1) = new_force(3*i-1) + (lj_force)*(y/r)
                new_force(3*i) = new_force(3*i) + (lj_force)*(z/r)
                new_force(3*j-2) = new_force(3*j-2) - (lj_force)*(x/r)
                new_force(3*j-1) = new_force(3*j-1) - (lj_force)*(y/r)
                new_force(3*j) = new_force(3*j) - (lj_force)*(z/r)
            end if 
        end do 
    end do 
    new_acc=new_force/dfloat(mass)

!UPDATE VELOCITY 

    ke= 0.0d0
    avr_vel_x = 0.0d0
    avr_vel_y = 0.0d0 
    avr_vel_z = 0.0d0 

    do i=1,n_part

        vel(3*i-2) = vel(3*i-2)+ (delta_t*0.5d0*(force(3*i-2)+ new_force(3*i-2)))
        vel(3*i-1) = vel(3*i-1)+ (delta_t*0.5d0*(force(3*i-1)+ new_force(3*i-1)))
        vel(3*i) = vel(3*i)+ (delta_t*0.5d0*(force(3*i)+ new_force(3*i)))

        avr_vel_x = avr_vel_x + vel(3*i-2)
        avr_vel_y = avr_vel_y + vel(3*i-1)
        avr_vel_z = avr_vel_z + vel(3*i)

        ke = ke + vel(3*i-2)*vel(3*i-2)
        ke = ke + vel(3*i-1)*vel(3*i-1)
        ke = ke + vel(3*i)*vel(3*i)
    
    end do 

    ke = 0.50d0*dfloat(mass)*ke 
    pot_energy = new_pot_energy 
    force=new_force 

!THERMOSTAT 
    if(mod(k,100)==0) then 
        theoryke = 1.5d0*dfloat(n_part)*kb_T
        scalef=dsqrt(theoryke/ke)
        vel=scalef*vel 

        ke=0.0d0 
        do i=1,n_part
            ke = ke + vel(3*i-2)*vel(3*i-2)
            ke = ke + vel(3*i-1)*vel(3*i-1)
            ke = ke + vel(3*i)*vel(3*i)
        end do 
        ke = 0.5*ke 
            write(*,*) 'ke',ke/dfloat(n_part)
            write(*,*) 'pe',(new_pot_energy)/dfloat(n_part)
            write(*,*) 'te',(ke+new_pot_energy)/dfloat(n_part)
    end if 

endsubroutine pos_update 

subroutine posinit
    use sim_para
    implicit none
    !Integer :: i,j
    Real*8 :: r1,r2,r3,distance
    Real*8 :: ran1,sq1,sq2,sq3,zz1
    zz1=ran1(zzz)
    !INITIALIZE THE POSITION
    pos(1)=dfloat(lx)*ran1(zzzz)
    pos(2)=dfloat(ly)*ran1(zzzz)
    pos(3)=dfloat(lz)*ran1(zzzz)
    do i=2,n_part 
11      pos(3*i-2)=dfloat(lx)*ran1(zzzz)
        pos(3*i-1)=dfloat(ly)*ran1(zzzz)
        pos(3*i)=dfloat(lz)*ran1(zzzz)
        do j=1,i-1
            r1=pos(3*i-2)-pos(3*j-2)
            r2=pos(3*i-1)-pos(3*j-1)
            r3=pos(3*i)-pos(3*j)

            r1=r1-nint(r1/dfloat(lx))*dfloat(lx)
            r2=r2-nint(r2/dfloat(ly))*dfloat(ly)
            r3=r3-nint(r3/dfloat(lz))*dfloat(lz)
            sq1=r1*r1
            sq2=r2*r2
            sq3=r3*r3
            distance=dsqrt(sq1+sq2+sq3)
            if (distance <= sigma) goto 11
        end do 
    end do 
    open(86,file='write_init_pos.dat',status='unknown')
    do i=1,n_part 
        write(86,*) pos(3*i-2),pos(3*i-1),pos(3*i)
    end do 

end subroutine posinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALCULATING INITIAL VELOCITY, FORCE AND ENERGY
subroutine vel_init 
    use sim_para
    implicit none

    real*8 :: vel_const,avg_vx,avg_vy,avg_vz 
    real*8 :: zz1,ran1 

    zz1 = ran1(zzz)
    vel_const=dsqrt(12.0d0)*dfloat(kb_T)

    !initialize the velocity
    do i=1,n_part
        vel(3*i-2)=vel_const*(ran1(zzzz)-0.5d0)
        vel(3*i-1)=vel_const*(ran1(zzzz)-0.5d0)
        vel(3*i)=vel_const*(ran1(zzzz)-0.5d0)
    end do 

    avg_vx=0.0d0; avg_vy=0.0d0; avg_vz=0.0d0 

    do i=1,n_part
        avg_vx = vel(3*i-2)+avg_vx; avg_vy = vel(3*i-1)+avg_vy 
        avg_vz = vel(3*i)+avg_vz 
    end do 

    avg_vx=avg_vx/dfloat(n_part)
    avg_vy=avg_vy/dfloat(n_part)
    avg_vz=avg_vz/dfloat(n_part)

    do i=1,n_part
        vel(3*i-2)=avg_vx-vel(3*i-2)
        vel(3*i-1)=avg_vy-vel(3*i-1)
        vel(3*i)=avg_vz-vel(3*i)
    end do 

    open(81,file='write_init_vel.dat',status='unknown')

    do i=1,n_part
        write(81,*) vel(3*i-2),vel(3*i-1),vel(3*i)
    end do 
endsubroutine vel_init 

Function ran1(idum)
    IMPLICIT NONE
    Integer,parameter :: K4B=selected_int_kind(9)
    Integer(K4B),Intent(INOUT) :: idum 
    Real*8 :: ran1 
    Integer(K4B),parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    Real,save :: am
    Integer(K4B),save :: ix=-1,iy=-1,k
    if(idum <= 0 .or. iy < 0) then
        am=nearest(1.0,-1.0)/IM 
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1
    end if 
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ 
    iy=IA*(iy-k*IQ)-IR*k 
    if (iy < 0) iy=iy+IM 
    ran1=am*ior(iand(IM,ieor(ix,iy)),1)
End Function ran1








