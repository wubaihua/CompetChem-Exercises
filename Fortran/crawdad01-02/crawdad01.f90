program crawdad01
    use LAPACK95
    
    implicit none
    
    integer ios,INFO
    integer:: natom
    integer ::i,j,k,l,n,lda,lwork,iwork,liwork
    real*8:: pi,xm,ym,zm,m,ixx,iyy,izz,ixy,ixz,iyz,tol
 !   real,allocatable :: x(:), y(:) ,z(:) ,nul(:)
!   real*8::mole(4,7)
 !   real*8::bond(7,7)
 !   real*8::ex(7,7),ey(7,7),ez(7,7),cosangel(7,7,7),angel(7,7,7),phi(7,7,7),theta(7,7,7,7)
    real*8::inertia(3,3),pmi(3),work(1)
   real*8,allocatable :: mole(:,:),bond(:,:),ex(:,:),ey(:,:),ez(:,:)
   real*8,allocatable :: cosangel(:,:,:),angel(:,:,:),phi(:,:,:),sintheta(:,:,:,:),theta(:,:,:,:)
    character*20:: name
    character*1::JOBZ,UPLO
   parameter(pi=3.1415926535897932)    
   
!***************************************************
!*******Read the molecular file:*****************
!***************************************************
    write(*,*)"please input the molecular :"
    read*, name
    
    open (10, file=trim(name)//'.txt', status='old')
    read(10,*) natom
    allocate(mole(4,natom))
    allocate(bond(natom,natom))
    allocate(ex(natom,natom))
    allocate(ey(natom,natom))
    allocate(ez(natom,natom))
    allocate(cosangel(natom,natom,natom))
    allocate(angel(natom,natom,natom))
    allocate(phi(natom,natom,natom))
    allocate(theta(natom,natom,natom,natom))
    allocate(sintheta(natom,natom,natom,natom))
    

    
    
    write(*,*)"The number of atoms:",natom
    do i=1,natom 
        read(10,*) mole(1:4,i)
        write(*,*) mole(1:4,i)
   end do
do
    read(10,iostat=ios)
     if (ios /= 0) exit
     end do
!write(*,*) A


!***************************************************
!**Calculate the bond lengths and bond matrix:******
!***************************************************   
write(*,*)"The bond matrix D(i,j):"
    do i=1,natom
        do j=1,natom
          bond(i,j)=sqrt((mole(2,i)-mole(2,j))**2+(mole(3,i)-mole(3,j))**2+(mole(4,i)-mole(4,j))**2)
        end do
     write(*,*)bond(i,1:7)
    end do
    write(*,*)"bond angel i-j-k:"
do i=1,natom
    do j=1,natom
        ex(i,j)=-(mole(2,i)-mole(2,j))/bond(i,j)
        ey(i,j)=-(mole(3,i)-mole(3,j))/bond(i,j)
        ez(i,j)=-(mole(4,i)-mole(4,j))/bond(i,j)
    end do
end do

!***************************************************
!**Calculate the bond angles:******
!***************************************************   

do i=3,natom
    do j=2,i-1
        do k=1,j-1
            if( bond(i,j)<4.0 .and. bond(k,j)<4.0 )then
    cosangel(i,j,k)=ex(j,i)*ex(j,k)+ey(j,i)*ey(j,k)+ez(j,i)*ez(j,k)
    phi(i,j,k)=acos(cosangel(i,j,k))
    angel(i,j,k)=(acos(cosangel(i,j,k)))*180/pi
   
    write(*,*)"angel",i,"-",j,"-",k,"=", angel(i,j,k)
    end if
        end do
    end do
end do
!***************************************************
!**Calculate the out-of-plane angles:******
!***************************************************   

write(*,*)"Out of plane angel i-j-k-l:"
do i=1,natom
    do k=1,natom
        do j=1,natom
            do l=1,j-1
                if(i/=j .and. i/=k .and. i/=l .and. k/=l .and. j/=k .and. bond(i,k)<4.0 .and. bond(k,j)<4.0 .and. bond(k,l)<4.0)then
                sintheta(i,j,k,l)=((ey(k,j)*ez(k,l)-ez(k,j)*ey(k,l))*ex(k,i)+(ez(k,j)*ex(k,l)-ex(k,j)*ez(k,l))*ey(k,i)+(ex(k,j)*ey(k,l)-ey(k,j)*ex(k,l))*ez(k,i))/sin(phi(j,k,l))
                if(sintheta(i,j,k,l)<-1.0)then
                    theta(i,j,k,l)=asin(-1.0)*180/acos(-1.0)
                else if(sintheta(i,j,k,l)>1.0)then
                    theta(i,j,k,l)=asin(1.0)*180/acos(-1.0)
                else
                    theta(i,j,k,l)=asin(sintheta(i,j,k,l))*180/acos(-1.0)
                end if
                
           write(*,*)"theta",i,"-",j,"-",k,"-",l,"=",theta(i,j,k,l)
                end if
                
            end do
        end do
    end do
end do
!***************************************************
!**Calculate the center of mass:********************
!***************************************************   
xm=0.0
ym=0.0
zm=0.0
m=0.0
do i=1,natom
    xm=xm+mole(1,i)*mole(2,i)
    ym=ym+mole(1,i)*mole(3,i)
    zm=zm+mole(1,i)*mole(4,i)
    m=m+mole(1,i)
end do
xm=xm/m
ym=ym/m
zm=zm/m
write(*,*)"The center of mass is ",xm,ym,zm
!***************************************************
!**Calculate the Moments of Inertia tensor:******
!***************************************************   
inertia(1,1)=0.0
inertia(2,2)=0.0
inertia(3,3)=0.0
inertia(1,2)=0.0
inertia(1,3)=0.0
inertia(2,3)=0.0
do i=1,natom
    inertia(1,1)=inertia(1,1)+mole(1,i)*(mole(3,i)**2+mole(4,i)**2)
    inertia(2,2)=inertia(2,2)+mole(1,i)*(mole(2,i)**2+mole(4,i)**2)
    inertia(3,3)=inertia(3,3)+mole(1,i)*(mole(2,i)**2+mole(3,i)**2)
    inertia(1,2)=inertia(1,2)+mole(1,i)*mole(2,i)*mole(3,i)
    inertia(2,3)=inertia(2,3)+mole(1,i)*mole(3,i)*mole(4,i)
    inertia(1,3)=inertia(1,3)+mole(1,i)*mole(2,i)*mole(4,i)
end do
!inertia(2,1)=inertia(1,2)
!inertia(3,2)=inertia(2,3)
!inertia(3,1)=inertia(1,3)
write(*,*) "The moment of inertia tensor is"
write(*,*) inertia(1:3,1)
write(*,*) inertia(1:3,2)
write(*,*) inertia(1:3,3)
!JOBZ='N'
!UPLO='U'
!n=3
!tol=1.0D-6

call dsyev('N','L',3,inertia,3,pmi,work,10,info)
!call dsyev(inertia,pmi)
write(*,*) "pmi="
write(*,*) pmi



 pause
end

Subroutine eigen_calc(a_op, a, a_matrix, n)

    Implicit None
    Integer n, ierr
    Real *8 a_op(n, n), a(n), a_matrix(n, n), number(3*n-1)

    !----------------------------------------------------------------------
    Call dsyev('V', 'L', n, a_op, n, a, number, 3*n-1, ierr) !Call the lapack to diagonalize
    !----------------------------------------------------------------------

    If (ierr/=0) Then !If it can't be Diagonalized, report error
        Write (*, *) 'Segmentation fault!'
    End If

    a_matrix(:, :) = a_op(:, :)

    End Subroutine eigen_calc
    !=======================================================================



!subroutine eigenvalue(A,n)
!implicit none
!integer,parameter::step=200
!integer::n,i,j,k
!real*8::A(n,n),x0(n),u(n),x1(n),lamda(n),au(n),err

!x0=1.0
!do i=1,n
!   u=x0/norm2(x0)
 !  call dot_matvec(A,u,n,x1)
 !  call dot_matvec(A,u,n,au)
 !  lamda(i)=dot_product(u,au)
 !  err=norm2(x1-x0)
!   if(err<1.0E-6)then
 !      exit

   
       


!end subroutine eigenvalue


!subroutine dot_matvec(A,x,n,b)
!integer::n,i,j
!real*8:: A(n,n),b(n),x(n)
!b=0.0
!do i=1,n
!    do j=1,n
!    b(i)=b(i)+A(i,j)*x(j)
!    end do
!end do
!end subroutine

!subroutine QR(n,A,Q,R)
!integer::n
!real*8::A(n,n),Q(n,n),R(n,n)
!
!
!
!
!
!
!
!
!
!
!
!
!end subroutine QR

!subroutine SOLVE(A,N,tezheng,tol)
!!QR分解特征值
!implicit real*8(a-z)
!integer::N
!real*8::A(N,N),tezheng(N)
!real*8::A1(N,N),Q(N,N),R(N,N)
!integer::i,j,k
!
!A1=A
!
!!循环迭代，最大允许迭代次数
!do i=1,200
!call GRAM_SC(A1,Q,R,N,N)
!A1=matmul(R,Q)
!!迭代停止标准
!do k=1,N
!ds=0d0
!ds=ds+A1(k,k)**2
!end do
!
!do j=1,N
!tezheng(j)=A1(j,j)
!end do
!
!if(ds < tol) exit
!end do
!end subroutine SOLVE
!
!subroutine GRAM_SC(A,Q,R,M,N)
!!采用修正的Gram-Schmidt分解进行QR分解
!implicit real*8(a-z)
!integer::M,N
!integer::i,j,k
!real*8::A(M,N),Q(M,N),R(N,N)
!real*8::vec_temp(M)
!
!R(1,1) = dsqrt(dot_product(A(:,1),A(:,1)))
!Q(:,1) = A(:,1)/R(1,1)
!do k=2,N
!do j=1,k-1
!R(j,k) = dot_product(Q(:,j),A(:,k))
!end do
!vec_temp = A(:,k)
!do j=1,k-1
!vec_temp = vec_temp - Q(:,j)*R(j,k)
!end do
!R(k,k)=dsqrt(dot_product(vec_temp,vec_temp))
!Q(:,k)=vec_temp/R(k,k)
!end do
!end subroutine GRAM_SC