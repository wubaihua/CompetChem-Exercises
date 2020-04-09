program crawdad04
    use LAPACK95
    implicit real*8(a-h,o-z)
    !real*8 s(7,7),sn(7,7),t(7,7),v(7,7),h(7,7),work(200),w(7),LS(7,7),Lambda(7,7),Sneg(7,7),F(7,7),eplsion(7),C(7,7),D(7,7),vt(7,7),Di(7,7),Cp(7,7)
    real*8,allocatable ::mole(:,:),s(:,:),sn(:,:),t(:,:),v(:,:),h(:,:),w(:),LS(:,:),work(:),miux(:,:),miuy(:,:),miuz(:,:)
    real*8,allocatable :: Lambda(:,:),Sneg(:,:),F(:,:),eplsion(:),C(:,:),D(:,:),vt(:,:),Di(:,:),Cp(:,:),eri(:,:,:,:),Fmo(:,:),erimo(:,:,:,:)
    !real*8 eri(7**2,7**2)
    !real*8 eri(7,7,7,7)
    !real*8 work(200)
    integer,external:: dindex
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
    allocate(erimo(nbasis,nbasis,nbasis,nbasis))
    
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
    call cpu_time(t1)
!***************************************************
!Compute the new Fock matrix and compute the new 
!density matrix and energies
!***************************************************      
 !F=H
500 write(*,*) "SCF in Step=",istep 
   do i=1,nbasis
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
    write(*,*) "The new Fock matrix is:"
    do i=1,nbasis
        write(*,*) F(i,:)
    end do
    
    do i=1,nbasis
        do j=1,i
            Cp(i,j)=F(i,j)
        end do
    end do
    call dsyev('V','L',nbasis,Cp,nbasis,eplsion,work,3*nbasis,info)   
    !do i=1,7
    !    eplsion(i,i)=w(i)
    !end do
    write(*,*)"The orbital energies are:"
    write(*,*) eplsion
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
 
!***************************************************
!Compute the error
!***************************************************     
    deltaE=Etoti-Etot
    rmsd=0.0
    do i=1,nbasis
        do j=1,nbasis
            rmsd=rmsd+(Di(i,j)-D(i,j))**2
        end do
    end do
    rmsd=rmsd**0.5
    write(*,*)"deltaE=",deltaE
    write(*,*)"RMSD=",rmsd
    if(deltaE<1.0E-6 .and. rmsd<1.0E-8)then
        write(*,*)"Hartree-Fock SCF Done!"
        call CPU_TIME(t2)
    else
        D=Di
        Etot=Etoti
        Eelec=Eeleci
        istep=istep+1
        goto 500
    end if
    write(*,*) "SCF steps=",istep
    write(*,*) "SCF time=",t2-t1
    
    
    D=Di
    Etot=Etoti
    Eelec=Eeleci
    write(*,*) "ESCF=",Etot
    
!***************************************************
!Transform the Two-Electron Integrals to the MO Basis
!Noddy Algorithm
!***************************************************  
   !write(*,*)"Noddy Algorithm"
   ! call cpu_time(t3)
   ! erimo=0.0
   ! do i=1,nbasis
   !     do j=1,nbasis
   !         do k=1,nbasis
   !             do l=1,nbasis
   !                 do i1=1,nbasis
   !                     do j1=1,nbasis
   !                         do k1=1,nbasis
   !                             do l1=1,nbasis
   !                                 erimo(i,j,k,l)=erimo(i,j,k,l)+C(i1,i)*C(j1,j)*C(k1,k)*C(l1,l)*eri(i1,j1,k1,l1)
   !                             end do
   !                         end do
   !                     end do
   !                 end do
   !             end do
   !         end do
   !     end do
   ! end do
    
  
!***************************************************
!Transform the Two-Electron Integrals to the MO Basis
!Smarter Algorithm 
!***************************************************      
    write(*,*)"Smarter Algorithm"
    call cpu_time(t3)
    erimo=0.0    
    !temp1=0.0
    !temp2=0.0
    !temp3=0.0
    !temp4=0.0
    sum1=0.0
    sum2=0.0
    sum3=0.0
    do i=1,nbasis
        do j=1,nbasis
            do k=1,nbasis
                do l=1,nbasis
                    do i1=1,nbasis
                        do j1=1,nbasis
                            do k1=1,nbasis
                                do l1=1,nbasis
                                    temp1=C(l1,l)*eri(i1,j1,k1,l1)
                                    sum1=temp1+sum1
                                end do
                                temp2=C(k1,k)*sum1
                                sum2=sum2+temp2
                                sum1=0.0
                            end do
                            temp3=C(j1,j)*sum2
                            sum3=sum3+temp3
                            sum2=0.0
                        end do
                        temp4=C(i1,i)*sum3
                        erimo(i,j,k,l)=erimo(i,j,k,l)+temp4
                        sum3=0.0
                    end do
                end do
            end do
        end do
    end do
    
                        
                    
                    
                    
                    
                    
!***************************************************
!Compute the MP2 energy
!***************************************************                          
    Emp2=0.0
    do i=1,nele/2
        do j=1,nele/2
            do ia=nele/2+1,nbasis
                do ib=nele/2+1,nbasis
                    Emp2=Emp2+erimo(i,ia,j,ib)*(2.0*erimo(i,ia,j,ib)-erimo(i,ib,j,ia))/(eplsion(i)+eplsion(j)-eplsion(ia)-eplsion(ib))
                end do
            end do
        end do
    end do
    write(*,*) "The MP2 energy is EMP2=",Emp2
    call cpu_time(t4)
    write(*,*) "MP2 time=",t4-t3
    write(*,*) "Etot=",Etot+Emp2
    
pause
end program
