program crawdad03
    use LAPACK95
    implicit real*8(a-h,o-z)
    !real*8 s(7,7),sn(7,7),t(7,7),v(7,7),h(7,7),work(200),w(7),LS(7,7),Lambda(7,7),Sneg(7,7),F(7,7),eplsion(7),C(7,7),D(7,7),vt(7,7),Di(7,7),Cp(7,7)
    real*8,allocatable ::mole(:,:),s(:,:),sn(:,:),t(:,:),v(:,:),h(:,:),w(:),LS(:,:),work(:),miux(:,:),miuy(:,:),miuz(:,:)
    real*8,allocatable :: Lambda(:,:),Sneg(:,:),F(:,:),eplsion(:),C(:,:),D(:,:),vt(:,:),Di(:,:),Cp(:,:),eri(:,:,:,:),Fmo(:,:)
    !real*8 eri(7**2,7**2)
    !real*8 eri(7,7,7,7)
    !real*8 work(200)
    real*8,allocatable ::Fi(:,:,:),ei(:,:,:),Dii(:,:,:)

    integer i,j,k,l,ij,ji,kl,lk,ijkl,jikl
 
!***************************************************
!Read the molecule geometry 
!***************************************************
   
    open (10,file="geom.txt",status='old')
    read(10,*) natom
    write(*,*) "natom=",natom
    allocate(mole(4,natom))
    write(*,*)"ZA----X---Y---Z"
    do i=1,natom 
        read(10,*) mole(1:4,i)
        write(*,*) mole(1:4,i)
    end do
    do
        read(10,iostat=ios)
        if (ios /= 0) exit
    end do
    
!***********************************************************
!Read the nbasis from s.txt and generate the matrices order
!***********************************************************
    open(19, file='s.txt', status='old')
    do 
        read(19,*,iostat=ierr) i,j,gin
        nbasis=i
        if(ierr<0)then
                exit
        end if
    end do
    write(*,*)"Nbasis=",nbasis
    nele=0
    do i=1,natom
        nele=nele+mole(1,i)
    end do
    write(*,*)"Nele=",nele
    
    allocate(s(nbasis,nbasis))
    allocate(sn(nbasis,nbasis))
    allocate(t(nbasis,nbasis))
    allocate(v(nbasis,nbasis))
    allocate(h(nbasis,nbasis))
    allocate(w(nbasis))
    allocate(LS(nbasis,nbasis))
    allocate(Lambda(nbasis,nbasis))
    allocate(Sneg(nbasis,nbasis))
    allocate(F(nbasis,nbasis))
    allocate(eplsion(nbasis))
    allocate(C(nbasis,nbasis))
    allocate(D(nbasis,nbasis))
    allocate(vt(nbasis,nbasis))
    allocate(Di(nbasis,nbasis))
    allocate(Cp(nbasis,nbasis))
    allocate(eri(nbasis,nbasis,nbasis,nbasis))
    allocate(work(3*nbasis))
    allocate(Fmo(nbasis,nbasis))
    allocate(miux(nbasis,nbasis))
    allocate(miuy(nbasis,nbasis))
    allocate(miuz(nbasis,nbasis))
    allocate(Fi(nbasis,nbasis,6))
    allocate(Dii(nbasis,nbasis,6))
    allocate(ei(nbasis,nbasis,6))
    
    
    
!***************************************************
!Read the nuclear repulsion energy=encu
!***************************************************
   
    open (11,file="enuc.txt",status='old')
    read(11,*) enuc
    write(*,*) "The nuclear repulsion energy=",enuc
!***************************************************
!Read the overlap integral matrix=s
!***************************************************  
    open(20, file='s.txt', status='old')
    do i=1,nbasis
        do j=1,i
        read(20,*) k,l,sn(i,j)
        end do
    end do
    s=sn
    do i=1,nbasis
        do j=i+1,nbasis
            s(i,j)=s(j,i)
        end do
    end do
    ! do i=1,7
    !    do j=1,i
    !    write(*,*) i,j,s(i,j)
    !    end do
    !end do
    write(*,*) "The overlap integral matrix is:"
    do i=1,nbasis
        write(*,*) s(:,i)
    end do
    
!***************************************************
!Read the kinetic integral matrix=t
!***************************************************  
    open(21, file='t.txt', status='old')
    do i=1,nbasis
        do j=1,i
        read(21,*) k,l,t(i,j)
        end do
    end do
    do i=1,nbasis
        do j=i+1,nbasis
            t(i,j)=t(j,i)
        end do
    end do
    write(*,*) "The kinetic integral matrix is:"
    do i=1,nbasis
        write(*,*) t(:,i)
    end do      
        
!***************************************************
!Read the nuclear-attraction integral matrix=v
!***************************************************  
    open(22, file='v.txt', status='old')
    do i=1,nbasis
        do j=1,i
        read(22,*) k,l,v(i,j)
        end do
    end do
    do i=1,nbasis
        do j=i+1,nbasis
            v(i,j)=v(j,i)
        end do
    end do
    write(*,*) "The nuclear-attraction integral matrix is:"
    do i=1,nbasis
        write(*,*) v(:,i)
   end do   
!***************************************************
!Calculate the core Hamiltonian=h
!***************************************************    
    h=t+v   
    write(*,*) "The core Hamiltonian matrix is:"
    do i=1,nbasis
        write(*,*) h(:,i)
   end do  
    
!***************************************************
!Read the two-electron repulsion integrals=eri
!***************************************************      
    eri=0.0
    open(23, file='eri.txt', status='old')
        do 
            read(23,*,iostat=ierr) i,j,k,l,din
            !ij=dindex(i,j)
            !kl=dindex(k,l)
            !ijkl=dindex(ij,kl)
            !ji=dindex(j,i)
            !lk=dindex(l,k)
            !jikl=ji*(ji+1)/2+kl  !            dindex(ji,kl)
            !ijlk=ij*(ij+1)/2+lk  !      dindex(ij,lk)
            !jilk=ji*(ji+1)/2+lk   !             dindex(ji,lk)
            !klij=kl*(kl+1)/2+ij    !            dindex(kl.ij)
            !klji=kl*(kl+1)/2+ji !dindex(kl.ji)
            !lkij=lk*(lk+1)/2+ij     !      dindex(lk.ij)
            !lkji=lk*(lk+1)/2+ji      !          dindex(lk.ji)
            !
            !eri(ijkl)=din
            !eri(ijlk)=eri(ijkl)
            !eri(jikl)=eri(ijkl)
            !eri(jilk)=eri(ijkl)
            !eri(klij)=eri(ijkl)
            !eri(klji)=eri(ijkl)
            !eri(lkij)=eri(ijkl)
            !eri(lkji)=eri(ijkl)
            
            eri(i,j,k,l)=din
            eri(i,j,l,k)=eri(i,j,k,l)
            eri(j,i,k,l)=eri(i,j,k,l)
            eri(j,i,l,k)=eri(i,j,k,l)
            eri(k,l,i,j)=eri(i,j,k,l)
            eri(k,l,j,i)=eri(i,j,k,l)
            eri(l,k,i,j)=eri(i,j,k,l)
            eri(l,k,j,i)=eri(i,j,k,l)
            
            
            
            !eri(dindex(i,j),dindex(k,l))=din
            !eri(dindex(j,i),dindex(k,l))=eri(dindex(i,j),dindex(k,l))
            !eri(dindex(i,j),dindex(l,k))=eri(dindex(i,j),dindex(k,l))
            !eri(dindex(j,i),dindex(l,k))=eri(dindex(i,j),dindex(k,l))
            !eri(dindex(k,l),dindex(i,j))=eri(dindex(i,j),dindex(k,l))
            !eri(dindex(k,l),dindex(j,i))=eri(dindex(i,j),dindex(k,l))
            !eri(dindex(l,k),dindex(i,j))=eri(dindex(i,j),dindex(k,l))
            !eri(dindex(l,k),dindex(j,i))=eri(dindex(i,j),dindex(k,l))
            
            
            
            if(ierr<0)then
                exit
            end if
        end do

        
        

        !write(*,*) "The two-electron repulsion integrals are:"
        !do i=1,7
        !    do j=1,7
        !        do k=1,7
        !            do l=1,7 
        !            write(*,*) i,j,k,l,eri(i,j,k,l)
        !            end do
        !        end do
        !    end do
        !end do
        !
        
        
!***************************************************
!Calculate S^-1/2=sneg
!***************************************************         
    call dsyev('V','L',nbasis,sn,nbasis,w,work,3*nbasis,info)   
    LS=sn 
    !write(*,*)"LS="
    !do i=1,nbasis
    !    write(*,*) LS(i,:)
    !end do
    !write(*,*) "Lambda="
    Lambda=0.0
    do i=1,nbasis
        Lambda(i,i)=w(i)**-0.5!Now Lambda is just ¦Ë^-1/2
        !write(*,*) Lambda(i,i)
    end do
    
    sneg=matmul(matmul(LS,Lambda),transpose(LS))
    write(*,*) "S^-1/2="  
        do i=1,nbasis
            write(*,*) sneg(:,i)
        end do
    
!***************************************************
!Build the initial density
!***************************************************         
    F=matmul(matmul(transpose(sneg),h),sneg)!Build the intital Fock matrix

    write(*,*) "The intital Fock matrix is:"
    do i=1,nbasis
        write(*,*) F(i,:)
    end do
    do i=1,nbasis
        do j=1,i
            C(i,j)=F(i,j)
        end do
    end do
    call dsyev('V','L',nbasis,C,nbasis,eplsion,work,3*nbasis,info)   
!call dgesvd('A', 'A',7, 7, F, 7, w, C, 10,vt, 10, work, 50, info)
    !do i=1,7
    !    eplsion(i,i)=w(i)
    !end do
    write(*,*)"The orbital energies are:"
    write(*,*) eplsion
    C=matmul(sneg,C)
    
    
    write(*,*) "The  initial MO coefficient matrix is:"
    do i=1,nbasis
        write(*,*) C(i,:)
    end do
    D=0.0
    do i=1,nbasis
        do j=1,nbasis
            do m=1,nele/2
                D(i,j)=D(i,j)+C(i,m)*C(j,m)
            end do
        end do
    end do
    write(*,*) "The  initial density matrix is:"
    do i=1,nbasis
        write(*,*) D(i,:)
    end do
    Dii(:,:,1)=D(:,:)
    
    !ei(:,:,1)=matmul(matmul(F,D),S)-matmul(matmul(S,D),F)
!***************************************************
!Calculate the SCF energy
!***************************************************        
    Eelec=0.0
    do i=1,nbasis
        do j=1,nbasis
            Eelec=Eelec+D(i,j)*(H(i,j)+F(i,j))
        end do
    end do
    
    
    
    
    !do i=1,5
    !    Eelec=Eelec+2*eplsion(i)
    !end do
    
    Etot=Eelec+enuc
    write(*,*)"The electronic energy is:",Eelec
    write(*,*)"The total energy is:",Etot

    
    istep=1
    Fi(:,:,1)=F(:,:)
    ei(:,:,istep)=matmul(matmul(F,D),S)-matmul(matmul(S,D),F)
    call cpu_time(t1)
!***************************************************
!Compute the new Fock matrix and compute the new 
!density matrix and energies
!***************************************************      
 !F=H
 write(*,*) "SCF in Step=",istep 
500   do i=1,nbasis
        do j=1,nbasis 
            F(i,j)=H(i,j)
            do k=1,nbasis
                do l=1,nbasis
                    !ij=i*(i+1)/2+j
                    !ij=dindex(i,j)
                    !kl=k*(k+1)/2+l
                    !ijkl=ij*(ij+1)/2+kl
                    !ik=i*(i+1)/2+k
                    !jl=j*(j+1)/2+l
                    !ikjl=ik*(ik+1)/2+jl
                    !F(i,j)=F(i,j)+D(k,l)*(2.0*eri(dindex(i,j),dindex(k,l))-eri(dindex(i,k),dindex(j,l)))
                    !F(i,j)=F(i,j)+D(k,l)*(2.0*eri(ijkl)-eri(ikjl))
                     F(i,j)=F(i,j)+D(k,l)*(2.0*eri(i,j,k,l)-eri(i,k,j,l))
                end do
            end do
        end do
    end do
    !write(*,*) "The new Fock matrix is:"
    !do i=1,nbasis
    !    write(*,*) F(i,:)
    !end do
    
530 do i=1,nbasis
        do j=1,i
            Cp(i,j)=F(i,j)
        end do
    end do
    call dsyev('V','L',nbasis,Cp,nbasis,eplsion,work,3*nbasis,info)   
    !do i=1,7
    !    eplsion(i,i)=w(i)
    !end do
    !write(*,*)"The orbital energies are:"
    !write(*,*) eplsion
    C=matmul(sneg,Cp)
    !write(*,*) "The  initial MO coefficient matrix is:"
    !do i=1,7
    !    write(*,*) C(i,:)
    !end do
    !do i=1,7
    !    do j=1,7
    !        Di(i,j)=0.0
    !    end do
    !end do
    Di=0.0
    do i=1,nbasis
        do j=1,nbasis
            do m=1,nele/2
                Di(i,j)=Di(i,j)+C(i,m)*C(j,m)
            end do
        end do
    end do
    !write(*,*) "The new density matrix is:"
    !do i=1,7
    !    write(*,*) Di(i,:)
    !end do
    Eeleci=0.0
    do i=1,nbasis
        do j=1,nbasis
            Eeleci=Eeleci+Di(j,i)*(H(i,j)+F(i,j))
        end do
    end do
    !do i=1,5
    !    Eeleci=Eeleci+2*eplsion(i)
    !end do
    Etoti=Eeleci+enuc
    write(*,*)"The electronic energy is:",Eeleci
    write(*,*)"The total energy is:",Etoti
    !F=0.0
    !if(istep>6)then
    !    do i=1,nbasis
    !        do j=1,nbasis 
    !            F(i,j)=H(i,j)
    !            do k=1,nbasis
    !                do l=1,nbasis
    !                    F(i,j)=F(i,j)+Di(k,l)*(2.0*eri(i,j,k,l)-eri(i,k,j,l))
    !                end do
    !            end do
    !        end do
    !    end do
    !end if
    
    
    istep=istep+1
    if(istep<=6)then
        ei(:,:,istep)=matmul(matmul(F,Di),S)-matmul(matmul(S,Di),F)
        em=0.0
        do i=1,nbasis
            do j=1,nbasis
                em=em+ei(i,j,istep)**2
            end do
        end do
        em=sqrt(em)
        write(*,*)"||ei||=",em
    else
        do i=1,5
            ei(:,:,i)=ei(:,:,i+1)
        end do
        ei(:,:,6)=matmul(matmul(F,Di),S)-matmul(matmul(S,Di),F)
        em=0.0
        do i=1,nbasis
            do j=1,nbasis
                em=em+ei(i,j,6)**2
            end do
        end do
        em=sqrt(em)
        write(*,*)"||ei||=",em
    end if
    
   
!***************************************************
!Compute the error
!***************************************************     
    deltaE=abs(Etoti-Etot)
    rmsd=0.0
    do i=1,nbasis
        do j=1,nbasis
            rmsd=rmsd+(Di(i,j)-D(i,j))**2
        end do
    end do 
    rmsd=rmsd**0.5
    
    write(*,*)"deltaE=",deltaE
    write(*,*)"RMSD=",rmsd
    if(deltaE<1.0E-8 .and. rmsd<1.0E-10)then
        write(*,*)"SCF Done!"
        call CPU_TIME(t2)
    else
        D=Di
        Etot=Etoti
        Eelec=Eeleci 
        if(istep>6)then
            do i=1,5
                Fi(:,:,i)=Fi(:,:,i+1)
            end do
            Fi(:,:,6)=F(:,:)
            !Dii(:,:,6)=D(:,:)
            write(*,*) "SCF steps=",istep
            write(*,*)"Using DIIS,order=6"
            call diis(nbasis,ei,Fi,F,6)
            goto 530
        else
            Fi(:,:,istep)=F(:,:)
            !Dii(:,:,istep)=D(:,:)
            write(*,*) "SCF steps=",istep
            !write(*,*)"Using DIIS,order=",istep
            !call diis(nbasis,ei,Fi,F,istep)
            goto 500
        end if
    end if
   
    write(*,*) "SCF time=",t2-t1
    
    D=Di
    Etot=Etoti
    Eelec=Eeleci
    

    
    
    
pause 

end program



subroutine diis(nbasis,ei,Fi,Fn,k)
    implicit real*8(a-h,o-z)
    real*8 Fi(nbasis,nbasis,k),Dii(nbasis,nbasis,k),ei(nbasis,nbasis,k),FS(nbasis,nbasis),Fn(nbasis,nbasis),D(nbasis,nbasis),e(nbasis,nbasis),S(nbasis,nbasis)
    real*8 B(k+1,k+1),c(k+1),bx(k+1)
    integer ipiv(k+1)
    B=0.0
    c=0.0
    !do istep=1,6
    !    F(:,:)=Fi(:,:,istep)
    !    D(:,:)=Dii(:,:,istep)
    !    ei(:,:,istep)=matmul(matmul(F,D),S)-matmul(matmul(S,D),F)
    !end do
    B(k+1,:)=-1.0
    B(:,k+1)=-1.0
    B(k+1,k+1)=0.0
    do i=1,k
        do j=1,k
            do m=1,nbasis
                do l=1,nbasis
                    B(i,j)=B(i,j)+ei(m,l,i)*ei(m,l,j)
                end do
            end do
        end do
    end do
    do i=1,k+1
        c(i)=0.0
        !do j=i+1,k+1
        !    B(i,j)=B(j,i)
        !end do
    end do
    c(k+1)=-1.0
    !write(*,*)"B=",B
    !write(*,*)"C=",c
    call dgesv(k+1,1,B,k+1,ipiv,c,k+1,info)
    Fn=0.0
    do i=1,k
        FS(:,:)=Fi(:,:,i)
        Fn=Fn+c(i)*FS
    end do
    !write(*,*)"Fn=",Fn
   
   
end subroutine