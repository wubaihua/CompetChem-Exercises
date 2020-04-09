program crawdad02
   
    use LAPACK95
   
    
    implicit none
    character*20::name
    integer::ios,natom,dim,i,j,k,l,step,m
    real*8,allocatable::mole(:,:),hessian(:,:),mwh(:,:),eighess(:),eigvec(:,:),C(:,:),work(:),z(:,:),wi(:)
    real*8::tol,scale,tau
    integer :: ilo,ihi,info
    real(kind=8) :: scalee(9)
    
    
    
    
    write(*,*)"please input the molecular:"
    read*, name
    open (10, file=trim(name)//'_geo.txt', status='old')
    read(10,*) natom
    allocate(mole(4,natom))
    write(*,*)"The number of atoms:",natom
    write(*,*)"The geo matrix is:"
    do i=1,natom 
        read(10,*) mole(1:4,i)
        write(*,*) mole(1:4,i)
    end do
    do
        read(10,iostat=ios)
        if (ios /= 0) exit
    end do   
    close(10)
    
    open (11, file=trim(name)//'_hessian.txt', status='old')
    read(11,*) natom
    dim=3*(natom**2)
    !allocate(hessian(3,dim))
    allocate(hessian(3*natom,3*natom))
    allocate(eighess(3*natom))
    allocate(eigvec(3*natom,3*natom))
    write(*,*)"The hessian matrix is:"
    do i=1,3*natom
        do j=1,natom
              read(11,*) hessian(i,3*j-2:3*j)
        end do
    end do
  
    do
        read(11,iostat=ios)
        if (ios /= 0) exit
    end do
    do i=1,3*natom
        write(*,*) hessian(i,1:3*natom),";" 
    end do
    
    
    
    
    
    allocate(mwh(3*natom,3*natom))
    do i=1,natom
        do k=3*i-2,3*i
            do j=1,natom
                do m=3*j-2,3*j
                    mwh(k,m)=hessian(k,m)/sqrt(mole(1,i)*mole(1,j))
                    mwh(k,m)=hessian(k,m)/sqrt(mole(1,i)*mole(1,j))
                    mwh(k,m)=hessian(k,m)/sqrt(mole(1,i)*mole(1,j))
                    !mwh(i,l)=hessian(j,m)/sqrt(mole(1,k)*mole(1,j))
            
                  
            end do
        end do
        end do
    end do
 write(*,*)"The Mass-Weight the Hessian Matrix is:"   
do l=1,3*natom
    write(*,*) mwh(l,1:3*natom),";"
end do





!allocate(C(3*natom,3*natom))
!call la_gemm(mwh,hessian,C)

!jobvl='N'
!jobvr='V'
!n=3*natom
!lda=2
!ldvl=1
!ldvr=2
!lwork=60
!wr=0.;wi=0.
eigvec=mwh

!write(*,*) "eig1"
!do i=1,3*natom
!write(*,*) eigvec(i,i)
!end do
allocate(work(3*natom))
allocate(z(3*natom,3*natom))
allocate(wi(3*natom))
call dgebal('N',3*natom,eigvec,3*natom,ilo,ihi,scale,info)
call dgehrd(3*natom,ilo,ihi,eigvec,3*natom,tau,work,-1,info)
call dhseqr('E','N',3*natom, ilo, ihi, eigvec, 3*natom,eighess,wi,z,1, work, 3*natom, info)
!call dgehrd(eigvec)
!call hseqr(eigvec,eighess)
write(*,*) "eig2:"
write(*,*) eighess

!call zgehrd(3*natom,1,3*natom,eigvec,3*natom,eighess,3*natom,ilo,ihi,scalee,info)
!write(*,*) "eig2"
!do i=1,3*natom
!write(*,*) eigvec(i,i)
!end do
!
!
!
!
!call SOLVE(mwh,3*natom,eighess,1.0D-6)
!write(*,*) "The eigenvalues of the mass-weighted Hessian is:"
!write(*,*) eighess
    
pause
end program 


subroutine SOLVE(A,N,tezheng,tol)
!QR分解特征值
implicit real*8(a-z)
integer::N
real*8::A(N,N),tezheng(N)
real*8::A1(N,N),Q(N,N),R(N,N)
integer::i,j,k

A1=A

!循环迭代，最大允许迭代次数
do i=1,200
call GRAM_SC(A1,Q,R,N,N)
A1=matmul(R,Q)
!迭代停止标准
do k=1,N
ds=0d0
ds=ds+A1(k,k)**2
end do

do j=1,N
tezheng(j)=A1(j,j)
end do

if(ds < tol) exit
end do
end subroutine SOLVE

subroutine GRAM_SC(A,Q,R,M,N)
!采用修正的Gram-Schmidt分解进行QR分解
implicit real*8(a-z)
integer::M,N
integer::i,j,k
real*8::A(M,N),Q(M,N),R(N,N)
real*8::vec_temp(M)

R(1,1) = dsqrt(dot_product(A(:,1),A(:,1)))
Q(:,1) = A(:,1)/R(1,1)
do k=2,N
do j=1,k-1
R(j,k) = dot_product(Q(:,j),A(:,k))
end do
vec_temp = A(:,k)
do j=1,k-1
vec_temp = vec_temp - Q(:,j)*R(j,k)
end do
R(k,k)=dsqrt(dot_product(vec_temp,vec_temp))
Q(:,k)=vec_temp/R(k,k)
end do
end subroutine GRAM_SC