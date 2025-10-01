program polar
use omp_lib
implicit none

integer lnum

doubleprecision, allocatable :: alm(:,:)
doubleprecision pi
integer l,m,i,j,mm,j1
doubleprecision k,randx
doubleprecision :: ascle(16388)
doubleprecision lpoint
doubleprecision lxxpoint,lxypoint,lyypoint,lypoint,lxpoint
doubleprecision lxxxpoint,lxxypoint,lxyypoint,lyyypoint
doubleprecision, allocatable :: lsum(:,:)
doubleprecision, allocatable :: lxsum(:,:), lysum(:,:)
doubleprecision, allocatable :: lxxsum(:,:), lxysum(:,:), lyysum(:,:)
doubleprecision, allocatable :: lxxxsum(:,:), lxxysum(:,:), lxyysum(:,:), lyyysum(:,:)
doubleprecision theta, phi
doubleprecision, allocatable :: func(:,:),q(:,:),u(:,:)
integer, allocatable :: dottyp(:,:)
doubleprecision det,disc
doubleprecision det00,disc00,det01,disc01,det10,disc10,det11,disc11
doubleprecision, allocatable :: q1(:,:),q2(:,:),u1(:,:),u2(:,:)
doubleprecision dpml,d2pml,d3pml
doubleprecision, allocatable :: spectrum(:)
doubleprecision saddle,beak,comet
doubleprecision msaddle,mbeak,mcomet
doubleprecision, allocatable :: sing(:,:)
integer singnum,sadnum,comnum,beanum
integer msadnum,mcomnum,mbeanum

doubleprecision q1p, q2p, u1p, u2p, reli, relj
doubleprecision timein,timeout,tim1, tim2
integer absi, absj,flag
doubleprecision norma1,norma2

doubleprecision, allocatable :: points1(:,:)

integer, allocatable :: mask(:,:)


character(len=10) :: mode, random
character(len=100) :: alm_filename,pointsname,areaname
integer :: sidenum, lmax,lmin
integer :: iostat
character(len=200) :: line
 open(unit=10, file='config.txt', status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        print *, 'file open error'
        stop
    end if
    
    ! Чтение и парсинг строк
    read(10, '(A)') line  ! Пропускаем разделитель
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) alm_filename
    
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) mode
    
   
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) random
    
    
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) sidenum
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) lmin
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) lmax
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) pointsname
    
    read(10, '(A)') line
    read(line(index(line, '=')+1:), *) areaname
    
    ! Закрытие файла
    close(10)
    
    ! Вывод результатов для проверки
    print *, 'Mode: ', trim(mode)
    print *, 'Gaussian or original?: ', trim(random)
    print *, 'Filename: ', trim(alm_filename)
    print *, 'Sidenum: ', sidenum
    print *, 'Lmax: ', lmax

if (lmin<2) then
 print*, 'l_max cannot be set to a value less than 2.'
 STOP
end if 

!call random_seed(put=[-778495513,2012632844, -1922164814,  -650205951 ,  705737051 ,-1089796552  , 512711102  ,   7839015&
!  ,-778495513,2012632844, -1922164814,  -650205951 ,  705737051 ,-1089796552  , 512711102  ,   7839015,&
!  -778495513,2012632844, -1922164814,  -650205951 ,  705757051 ,-1089796552  , 512711102  ,   7839015&
!  ,-778495513,2012632844, -1922164814,  -650205951 ,  705637051 ,-1089796552  , 512711102  ,   7839015,&
!  -778495513,2012632844, -1922164814,  -650205951 ,  705787051 ,-1089796552  , 512711102  ,   7839015&
!  ,-778495513,2012632844, -1922164814,  -650205951 ,  705837051 ,-1089796552  , 512711102  ,   7839015])
tim1=omp_get_wtime()


pi=4.d0*datan(1.d0)





lnum=sidenum
allocate(points1(3,5000000))
allocate(spectrum(lnum))

allocate(alm(lnum,lnum*2+1))


allocate(lxsum(2*lnum+1,lnum-1))
allocate(lysum(2*lnum+1,lnum-1))
allocate(lxxsum(2*lnum+1,lnum-1))
allocate(lxysum(2*lnum+1,lnum-1))
allocate(lyysum(2*lnum+1,lnum-1))

allocate(lxxxsum(2*lnum+1,lnum-1))
allocate(lxxysum(2*lnum+1,lnum-1))
allocate(lxyysum(2*lnum+1,lnum-1))
allocate(lyyysum(2*lnum+1,lnum-1))

allocate(dottyp(lnum-1,lnum*2))
allocate(q(lnum-1,lnum*2))
allocate(u(lnum-1,lnum*2))
allocate(q1(lnum-1,lnum*2))
allocate(u1(lnum-1,lnum*2))
allocate(q2(lnum-1,lnum*2))
allocate(u2(lnum-1,lnum*2))






allocate(mask(lnum-1,lnum*2))
open(1,file=trim(alm_filename),form="formatted",ACCESS='SEQUENTIAL', status='old'&
    ,action='read')
    
    do l=0,2048
        do m=-l,l
            read(1,*)k

            if (l>0) then
            if (m==0) then
                alm(l,m+lnum+1)=(-1)**m*k*dble(10)**8
                
            else
                alm(l,m+lnum+1)=(-1)**m*k*sqrt(2d0)*dble(10)**8

            end if
            end if
            
           
        end do
    end do
close(1)



    
spectrum=0d0
do l=2,lnum
    do m=-l,l
        spectrum(l)=spectrum(l)+alm(l,m+lnum+1)**2/(2*l+1.d0)
    end do
end do

if (trim(random) == 'gaussian'  ) then
do l=2,lnum
    lpoint=0
    do m=-l,l
        
            call random_std(randx)
            alm(l,m+lnum+1)=randx
            lpoint=lpoint+randx**2


    end do
    

     lypoint=dsqrt(spectrum(l)*(2*l+1.d0)/(lpoint))

     alm(l,:)=alm(l,:)*lypoint


end do
 end if

mask=1



print*,lnum
j=0

q=0
u=0
q1=0
q2=0
u1=0
u2=0



do l=2, lmax
norma1=l**2/dsqrt((l-1d0)*l*(l+1d0)*(l+2d0))
norma2=dble(l)**3/dsqrt((l-1d0)*l*(l+1d0)*(l+2d0))/2048d0



lxsum=0
lysum=0
lxxsum=0
lxysum=0
lyysum=0
lxxxsum=0
lxxysum=0
lxyysum=0
lyyysum=0
 

timein=omp_get_wtime()
!$omp parallel
!$omp do private (i,m,theta,ascle,lpoint,lxpoint,lypoint,lxxpoint,lxypoint,lyypoint,lxyypoint,lxxxpoint,lxxypoint,lyyypoint)&
!$omp & private (dpml,d2pml,d3pml,mm)&
!$omp & private (det,disc,q1p,q2p,u1p,u2p)&
!$omp & schedule (dynamic) 
do i=2, lnum
    
    
    theta=(i-1d0)*pi/(lnum)
    
    call asc_legtrid(l,theta,ascle)
    

    do m=-l,-1
        
  

           
        mm=abs(m)
        
        dpml=dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)*ascle(mm+1+1)+mm/&
        dtan(theta)&
        /dble(l)*ascle(mm+1)
        
                    
        d2pml= (-dble(l-mm)*dble(l+mm+1) + dble(mm)**2/dtan(theta)**2 - dble(mm)/dsin(theta)**2)/dble(l)**2 * ascle(mm+1) &
        - dsqrt(dble(l-mm)*dble(l+mm+1))/dtan(theta)/dble(l)**2 * ascle(mm+2)
        

        
                            
        d3pml=(dsqrt(dble(l-mm)*dble(l-mm-1)*dble(l-mm-2)*dble(l+mm+1)*dble(l+mm+2)*dble(l+mm+3)))&
        /dble(l)**3*ascle((mm+3)+1)+&
        (3d0*dble(mm)+3)/tan(theta)*dsqrt(dble(l-mm)*dble(l-mm-1)*dble(l+mm+1)*dble(l+mm+2))&
        /dble(l)**3*ascle((mm+2)+1)+&
        ((3*dble(mm)**2+3*dble(mm)+1)/(tan(theta))**2-(3*dble(mm)+1)/(sin(theta))**2)*dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)**3*ascle((mm+1)+1)+&
        (dble(mm)**3/(tan(theta))**3-dble(mm)**2/tan(theta)/(sin(theta))**2+&
        2*dble(mm)/(tan(theta)*sin(theta)**2)*(1d0-mm))&
        /dble(l)**3*ascle((mm)+1)
        
        
        
        

        
        d3pml=-(dble(l-mm-1)*dble(l+mm+2)-(dble(mm)**2+3*mm+3)&
        /dtan(theta)**2+dble(3*mm+1)/sin(theta)**2)&
        *dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)**3*ascle((mm+1)+1)+&
        ((2*mm-3*dble(mm)**2)/sin(theta)**2-dble(mm-1)*dble(l-mm)*dble(l+mm+1)+dble(mm)**3/dtan(theta)**2)&
        /dtan(theta)/dble(l)**3*ascle((mm)+1)
                

                
                                   
            lxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml!sin
            lypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*ascle((mm)+1)*dble(mm)/dble(l)!cos
            lxxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d2pml!sin
            lyypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*ascle((mm)+1)*(mm/dble(l))**2*(-1)!sin
            lxypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml*abs(m/dble(l))!cos
                    
            lxxxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d3pml!sin
            lxxypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d2pml*(mm/dble(l))!cos
            lxyypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml*(mm/dble(l))**2*(-1)!sin
            lyyypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*ascle(mm+1)*(-(mm/dble(l))**3)!cos
                    
          !print*,d3pml,lxxxpoint,m          
                    
            lxsum(2*(-m)+2,i-1)= lxsum(2*(-m)+2,i-1)+lxpoint/sqrt(2.d0)
            
            lysum(2*(-m)+1,i-1)= lysum(2*(-m)+1,i-1)+lypoint/sqrt(2.d0)
            
            lxxsum(2*(-m)+2,i-1)= lxxsum(2*(-m)+2,i-1)+lxxpoint/sqrt(2.d0)
            
            lxysum(2*(-m)+1,i-1)= lxysum(2*(-m)+1,i-1)+lxypoint/sqrt(2.d0)
            
            lyysum(2*(-m)+2,i-1)= lyysum(2*(-m)+2,i-1)+lyypoint/sqrt(2.d0)
            
            lxxxsum(2*(-m)+2,i-1)= lxxxsum(2*(-m)+2,i-1)+lxxxpoint/sqrt(2.d0)
            
            lxxysum(2*(-m)+1,i-1)= lxxysum(2*(-m)+1,i-1)+lxxypoint/sqrt(2.d0)
            
            lxyysum(2*(-m)+2,i-1)= lxyysum(2*(-m)+2,i-1)+lxyypoint/sqrt(2.d0)
            
            lyyysum(2*(-m)+1,i-1)= lyyysum(2*(-m)+1,i-1)+lyyypoint/sqrt(2.d0)
      end do
      
      mm=0
      m=0
        
        dpml=dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)*ascle(mm+1+1)+mm/&
        dtan(theta)&
        /dble(l)*ascle(mm+1)
        
                    
                    
        d2pml = (-dble(l-mm)*dble(l+mm+1) + dble(mm)**2/dtan(theta)**2 - dble(mm)/dsin(theta)**2)/dble(l)**2 * ascle(mm+1) &
        - dsqrt(dble(l-mm)*dble(l+mm+1))/dtan(theta)/dble(l)**2 * ascle(mm+2)
        
                            
        d3pml=-(dble(l-mm-1)*dble(l+mm+2)-(dble(mm)**2+3*mm+3)&
        /dtan(theta)**2+dble(3*mm+1)/sin(theta)**2)&
        *dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)**3*ascle((mm+1)+1)+&
        ((2*mm-3*dble(mm)**2)/sin(theta)**2-dble(mm-1)*dble(l-mm)*dble(l+mm+1)+dble(mm)**3/dtan(theta)**2)&
        /dtan(theta)/dble(l)**3*ascle((mm)+1)
        
      
            lxxxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d3pml
            lxxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d2pml
            lxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml
                    
            lxxxsum(1,i-1)=lxxxsum(1,i-1)+lxxxpoint         
            lxxsum(1,i-1)=lxxsum(1,i-1)+lxxpoint
            lxsum(1,i-1)=lxsum(1,i-1)+lxpoint
                    
        do m=1,l
        
            mm=m
        
        dpml=dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)*ascle(mm+1+1)+mm/&
        dtan(theta)&
        /dble(l)*ascle(mm+1)
        
                                       
        d2pml = (-dble(l-mm)*dble(l+mm+1) + dble(mm)**2/dtan(theta)**2 - dble(mm)/dsin(theta)**2)/dble(l)**2 * ascle(mm+1) &
        - dsqrt(dble(l-mm)*dble(l+mm+1))/dtan(theta)/dble(l)**2 * ascle(mm+2)
        
                            
        d3pml=-(dble(l-mm-1)*dble(l+mm+2)-(dble(mm)**2+3*mm+3)&
        /dtan(theta)**2+dble(3*mm+1)/sin(theta)**2)&
        *dsqrt(dble(l-mm)*dble(l+mm+1))&
        /dble(l)**3*ascle((mm+1)+1)+&
        ((2*mm-3*dble(mm)**2)/sin(theta)**2-dble(mm-1)*dble(l-mm)*dble(l+mm+1)+dble(mm)**3/dtan(theta)**2)&
        /dtan(theta)/dble(l)**3*ascle((mm)+1)
             
        
                    lxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml!*cos(m*phi(j))
                    lypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*ascle(mm+1)*(-m/dble(l))!*sin(m*phi(j))
                    lxxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d2pml!*cos(m*phi(j))
                    lyypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*ascle(m+1)*(mm/dble(l))**2*(-1)!cos(m*phi(j)))
                    lxypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml*(-mm/dble(l))!*sin(m*phi(j))
                    
                    lxxxpoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d3pml!cos
                    lxxypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*d2pml*(-mm/dble(l))!sin
                    lxyypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*dpml*(mm/dble(l))**2*(-1)!cos
                    lyyypoint=alm(l,m+lnum+1)*dsqrt(1d0/2d0/pi)*ascle(mm+1)*((mm/dble(l))**3)!sin
                    
            lxsum(2*(m)+1,i-1)=lxsum(2*(m)+1,i-1)+lxpoint/sqrt(2.d0)
            
            lysum(2*(m)+2,i-1)=lysum(2*(m)+2,i-1)+lypoint/sqrt(2.d0)
            
            lxxsum(2*(m)+1,i-1)=lxxsum(2*(m)+1,i-1)+lxxpoint/sqrt(2.d0)
           
            lxysum(2*(m)+2,i-1)=lxysum(2*(m)+2,i-1)+lxypoint/sqrt(2.d0)
           
            lyysum(2*(m)+1,i-1)=lyysum(2*(m)+1,i-1)+lyypoint/sqrt(2.d0)
            
            
            lxxxsum(2*(m)+1,i-1)=lxxxsum(2*(m)+1,i-1)+lxxxpoint/sqrt(2.d0)
            
            lxxysum(2*(m)+2,i-1)=lxxysum(2*(m)+2,i-1)+lxxypoint/sqrt(2.d0)
            
            lxyysum(2*(m)+1,i-1)=lxyysum(2*(m)+1,i-1)+lxyypoint/sqrt(2.d0)
            
            lyyysum(2*(m)+2,i-1)=lyyysum(2*(m)+2,i-1)+lyyypoint/sqrt(2.d0)
       
         
       
    end do
    
    lxsum(2,i-1)=lxsum(2*(lnum)+1,i-1)
    lysum(2,i-1)=lysum(2*(lnum)+1,i-1)
    lxxsum(2,i-1)=lxxsum(2*(lnum)+1,i-1)
    lxysum(2,i-1)=lxysum(2*(lnum)+1,i-1)
    lyysum(2,i-1)=lyysum(2*(lnum)+1,i-1)
    lxxysum(2,i-1)=lxxysum(2*(lnum)+1,i-1)
    lxxxsum(2,i-1)=lxxxsum(2*(lnum)+1,i-1)
    lxyysum(2,i-1)=lxyysum(2*(lnum)+1,i-1)
    lyyysum(2,i-1)=lyyysum(2*(lnum)+1,i-1)




    theta=(i-1d0)*pi/(lnum)
        call realft(lxsum(:2*lnum,i-1),2*lnum,-1)
        call realft(lysum(:2*lnum,i-1),2*lnum,-1)
        call realft(lxxsum(:2*lnum,i-1),2*lnum,-1)
        call realft(lxysum(:2*lnum,i-1),2*lnum,-1)
        call realft(lyysum(:2*lnum,i-1),2*lnum,-1)
        call realft(lyyysum(:2*lnum,i-1),2*lnum,-1)
        call realft(lxxxsum(:2*lnum,i-1),2*lnum,-1)
        call realft(lxxysum(:2*lnum,i-1),2*lnum,-1)
        call realft(lxyysum(:2*lnum,i-1),2*lnum,-1)

        do j=1,lnum*2




        if (trim(mode) == 'E') then    
            q(i-1,j)=q(i-1,j)+(lxxsum(j,i-1)-1d0/tan(theta)/dble(l)*lxsum(j,i-1)-&
            1.d0/(sin(theta))**2*lyysum(j,i-1))*norma1
            u(i-1,j)=u(i-1,j)+2.0/sin(theta)*(&
            lxysum(j,i-1)-1d0/tan(theta)/dble(l)*lysum(j,i-1))*norma1!2.0/sin(theta(i))*(fxypoint)
            
            q1(i-1,j)=q1(i-1,j)+(lxxxsum(j,i-1)&
            -1.d0/((sin(theta))**2)*lxyysum(j,i-1)&
            -1d0/tan(theta)/dble(l)*lxxsum(j,i-1)&
            +1d0/((sin(theta))**2)/dble(l)**2*lxsum(j,i-1)&
            +2d0/tan(theta)/((sin(theta))**2)/dble(l)*lyysum(j,i-1))*norma2
            
            u1(i-1,j)=u1(i-1,j)+&
            2d0/sin(theta)*(lxxysum(j,i-1)&
            -2d0/tan(theta)/dble(l)*lxysum(j,i-1)&
            +((1d0/tan(theta))**2+1.d0/(sin(theta))**2)/dble(l)**2*lysum(j,i-1))&
            *norma2
            
            q2(i-1,j)=q2(i-1,j)+(lxxysum(j,i-1)/sin(theta)&
            -5d0/tan(theta)/sin(theta)/dble(l)*lxysum(j,i-1)&
            +4d0/(tan(theta))**2/sin(theta)/dble(l)**2*lysum(j,i-1)&
            -lyyysum(j,i-1)/(sin(theta))**3)*norma2
            
            u2(i-1,j)=u2(i-1,j)+2*(lxyysum(j,i-1)/(sin(theta))**2&
            +1d0/tan(theta)/dble(l)*lxxsum(j,i-1)&
            -2d0/tan(theta)/(sin(theta))**2/dble(l)*lyysum(j,i-1)&
            -1d0/(tan(theta))**2/dble(l)**2*lxsum(j,i-1))*norma2
        else
            u(i-1,j)=u(i-1,j)+(lxxsum(j,i-1)-1d0/tan(theta)/dble(l)*lxsum(j,i-1)-&
            1.d0/(sin(theta))**2*lyysum(j,i-1))*norma1
            q(i-1,j)=q(i-1,j)-2.0/sin(theta)*(&
            lxysum(j,i-1)-1d0/tan(theta)/dble(l)*lysum(j,i-1))*norma1!2.0/sin(theta(i))*(fxypoint)
            
            u1(i-1,j)=u1(i-1,j)+(lxxxsum(j,i-1)&
            -1.d0/((sin(theta))**2)*lxyysum(j,i-1)&
            -1d0/tan(theta)/dble(l)*lxxsum(j,i-1)&
            +1d0/((sin(theta))**2)/dble(l)**2*lxsum(j,i-1)&
            +2d0/tan(theta)/((sin(theta))**2)/dble(l)*lyysum(j,i-1))*norma2
            
            q1(i-1,j)=q1(i-1,j)-&
            2d0/sin(theta)*(lxxysum(j,i-1)&
            -2d0/tan(theta)/dble(l)*lxysum(j,i-1)&
            +((1d0/tan(theta))**2+1.d0/(sin(theta))**2)/dble(l)**2*lysum(j,i-1))&
            *norma2
            
            u2(i-1,j)=u2(i-1,j)+(lxxysum(j,i-1)/sin(theta)&
            -5d0/tan(theta)/sin(theta)/dble(l)*lxysum(j,i-1)&
            +4d0/(tan(theta))**2/sin(theta)/dble(l)**2*lysum(j,i-1)&
            -lyyysum(j,i-1)/(sin(theta))**3)*norma2
            
            q2(i-1,j)=q2(i-1,j)-2*(lxyysum(j,i-1)/(sin(theta))**2&
            +1d0/tan(theta)/dble(l)*lxxsum(j,i-1)&
            -2d0/tan(theta)/(sin(theta))**2/dble(l)*lyysum(j,i-1)&
            -1d0/(tan(theta))**2/dble(l)**2*lxsum(j,i-1))*norma2
        
        end if
        
        q1p=q1(i-1,j)
            q2p=q2(i-1,j)
            u1p=u1(i-1,j)
            u2p=u2(i-1,j)
        det=q1p*u2p-q2p*u1p
        disc=4*(u1p+2*q2p)**3*u1p+(u1p+2*q2p)**2*(2*q1p-u2p)**2&
        -4*u2p*(2*q1p-u2p)**3-18*u2p*(u1p+2*q2p)*&
        (2*q1p-u2p)*u1p-27*u2p**2*u1p**2
        
        if (det>0) then
                if (disc>0) then!beak
                    dottyp(i-1,j)=3
                else!comet
                    dottyp(i-1,j)=2
                end if
        end if
        if (det<0) then
                if (disc>0) then!saddle
                    dottyp(i-1,j)=1
                else
                    dottyp(i-1,j)=40
                end if
        end if

            
        end do
end do

!$omp end do
!$omp end parallel
timeout=omp_get_wtime()
print*,l, 'l',timeout-timein,'time'

singnum=0
do j=1,lnum*2
 do i=2,lnum
  theta=(i-1d0)*pi/(lnum)
  j1=mod(j,lnum*2)+1   

  call findsingu(j,i,lnum, q(i-1,j),q(i-1,j1),q(i,j),q(i,j1),&
  u(i-1,j),u(i-1,j1),u(i,j),u(i,j1),&
  q2(i-1,j),q2(i-1,j1),q2(i,j),q2(i,j1),&
  u2(i-1,j),u2(i-1,j1),u2(i,j),u2(i,j1),&
  q1(i-1,j),q1(i-1,j1),q1(i,j),q1(i,j1),&
  u1(i-1,j),u1(i-1,j1),u1(i,j),u1(i,j1),&
  q2p,q1p,u2p,u1p,reli,relj,flag&
  )
  det=q1p*u2p-q2p*u1p
  disc=4*(u1p+2*q2p)**3*u1p+(u1p+2*q2p)**2*(2*q1p-u2p)**2&
  -4*u2p*(2*q1p-u2p)**3-18*u2p*(u1p+2*q2p)*&
  (2*q1p-u2p)*u1p-27*u2p**2*u1p**2
  if (flag==1) then
   singnum=singnum+1
   points1(1,singnum)=(reli-1)*pi/(lnum)
   points1(2,singnum)=(relj-1)*pi/(lnum)
   if (det>0) then
    if (disc>0) then!beak
    
       points1(3,singnum)=2

    else!comet

       points1(3,singnum)=3
    end if
   end if
    if (det<0) then
     if (disc>0) then!saddle
    

       points1(3,singnum)=1
     else

     end if
   end if
  end if
 end do
end do



if (l==lmax) then
open (unit=10,form="formatted",ACCESS='SEQUENTIAL',file=trim(pointsname),action="write")
open (unit=55,form="formatted",ACCESS='SEQUENTIAL',file=trim(areaname),action="write")

do j=1,lnum*2
    write (55,*) dottyp(:,j)
end do
do i=1,singnum
 write (10,*) points1(:,i)
end do
close(10)
close(55)
end if

enddo

tim2=omp_get_wtime()

print*,'computation time',tim2-tim1

contains 




function ind(l_,m_)
    implicit none
    integer, intent(in) :: l_,m_
    integer ind
    ind=(l_+1)*l_/2+m_+1
    return
end function




SUBROUTINE realft(data,n,isign)
    implicit none
INTEGER isign,n
DOUBLE PRECISION data(n)


INTEGER i,i1,i2,i3,i4,n2p3
double precision c1,c2,h1i,h1r,h2i,h2r,wis,wrs
DOUBLE PRECISION theta,wi,wpi,wpr, wr,wtemp 
theta=4.d0*datan(1.d0)/dble(n/2) 
c1=0.5
if (isign.eq.1) then
  c2=-0.5
  call four1(data,n/2,+1) 
else
  c2=0.5 
  theta=-theta
endif
wpr=-2.0d0*sin(0.5d0*theta)**2
wpi=sin(theta)
wr=1.0d0+wpr
wi=wpi
n2p3=n+3
do i=2,n/4 
  i1=2*i-1
  i2=i1+1
  i3=n2p3-i2
  i4=i3+1
  wrs=dble(wr)
  wis=dble(wi)
  h1r=c1*(data(i1)+data(i3)) 
  h1i=c1*(data(i2)-data(i4))
  h2r=-c2*(data(i2)+data(i4))
  h2i=c2*(data(i1)-data(i3))
  data(i1)=h1r+wrs*h2r-wis*h2i 
  data(i2)=h1i+wrs*h2i+wis*h2r
  data(i3)=h1r-wrs*h2r+wis*h2i
  data(i4)=-h1i+wrs*h2i+wis*h2r
  wtemp=wr
  wr=wr*wpr-wi*wpi+wr
  wi=wi*wpr+wtemp*wpi+wi
enddo 

if (isign.eq.1) then
  h1r=data(1)
  data(1)=h1r+data(2)
  data(2)=h1r-data(2) 
else
  h1r=data(1)
  data(1)=c1*(h1r+data(2))
  data(2)=c1*(h1r-data(2))
  call four1(data,n/2,-1) 
endif
return
END


SUBROUTINE four1(dat,nn,isign1)
    implicit none
INTEGER isign1,nn
doubleprecision dat(2*nn)
!Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
!data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as −1.
!data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
!MUST be an integer power of 2 (this is not checked for!).
INTEGER i,istep,j,m,mmax,n
doubleprecision tempi,tempr
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp, pi !Double precision for the trigonometric recurrences.
pi=4.d0*datan(1.d0)

n=2*nn
j=1

do  i=1,n,2 !This is the bit-reversal section of the routine.
    if(j>i)then
        tempr=dat(j) !Exchange the two complex numbers.
        tempi=dat(j+1)
        dat(j)=dat(i)
        dat(j+1)=dat(i+1)
        dat(i)=tempr
        dat(i+1)=tempi
    endif
    m=nn
    do while ((m>=2).and.(j>m))
        j=j-m
        m=m/2
    end do
    j=j+m
enddo
mmax=2 !Here begins the Danielson-Lanczos section of the routine.
do while (n>mmax) !Outer loop executed log2 nn times.
    istep=2*mmax
    theta=2.d0*pi/(isign1*mmax) !Initialize for the trigonometric recurrence.
    wpr=-2.d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.d0
    wi=0.d0
    do m=1,mmax,2 !Here are the two nested inner loops

        do i=m,n,istep
            j=i+mmax !This is the Danielson-Lanczos formula:
            tempr=wr*dat(j)-wi*dat(j+1)
            tempi=wr*dat(j+1)+wi*dat(j)
            dat(j)=dat(i)-tempr
            dat(j+1)=dat(i+1)-tempi
            dat(i)=dat(i)+tempr
            dat(i+1)=dat(i+1)+tempi
        enddo
        wtemp=wr !Trigonometric recurrence.
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
    enddo
    mmax=istep
end do !All done.
end

subroutine random_stduniform(u)
   implicit none
   doubleprecision,intent(out) :: u
   doubleprecision :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform
subroutine random_std(x)
   implicit none
   doubleprecision,intent(out) :: x
   doubleprecision :: pi
   doubleprecision :: u1,u2
   pi=4.d0*datan(1.d0)
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_std

subroutine asc_legtrid(LL_,theta_,p_)!NN<=256 LL<=1030
    implicit none
    doubleprecision pi_
    doubleprecision, intent(in) :: theta_
    integer, intent(in) :: LL_
    integer m_,l_

    doubleprecision, intent(out) :: p_(1:16388)
    doubleprecision :: pl0_(1:16385),b(16383),c(16383),f(16383)
    doubleprecision x_,y_
    doubleprecision coef_, c1_,c2_, coef1_, coef2_
    pi_=4.d0*datan(1.d0)

    !LL=32
    !p1_=0.d0
    !вводим новое представление полиномов Лежандра
    !было: p(m+1,l+1,j)
    !будет: если m>l то 0, инчае p((l+1)*l/2+m+1,j)
    p_=0
    x_=cos(theta_)
    y_=sin(theta_)
    
    p_(0+1)=1.d0/dsqrt(2.d0)
    pl0_(0+1)=1.d0/dsqrt(2.d0)
    if ((theta_.ne.pi_/2d0).and. l .ne. 2) then
    do m_=1,LL_
        coef_=-dsqrt((2.d0*dble(m_)+1.d0)/(2.d0*dble(m_)))

        p_(m_+1)=coef_*p_(m_-1+1)*y_
            !print*,x(j)

    end do
    
    m_=1
    coef_=sqrt(2d0*m_)
    pl0_(1+1)=-coef_*x_/y_*p_(m_+1)
    
    
    m_=0
    do  l_=2,LL_

            c1_=(2.d0*dble(l_)+1.d0)*(2.d0*dble(l_)-1.d0)&
                 /(dble(l_)+dble(m_))/(dble(l_)-dble(m_))
            coef1_=dsqrt(c1_)

            c2_=(dble(l_)-1.d0-dble(m_))*(dble(l_)-1.d0+dble(m_))*&
              (2.d0*dble(l_)+1.d0)&
              /(dble(l_)+dble(m_))/(dble(l_)-dble(m_))/(2.d0*dble(l_)-3.d0)
            coef2_=dsqrt(c2_)


            pl0_(l_+1)=coef1_*x_*pl0_(l_-1+1)-coef2_*pl0_(l_-2+1)

    end do

    l_=LL_
    
    
        m_=0

            c1_=4d0*(dble(m_)+1.d0)**2&
                 /dble(l_-m_)/(l_+m_+1d0)
            coef1_=dsqrt(c1_)

            c2_=(l_+m_+2d0)*(l_-m_-1d0)&
              /dble(l_-m_)/(l_+m_+1d0)
            coef2_=dsqrt(c2_)

            
            b(m_+1)=coef1_*x_/y_
            c(m_+1)=coef2_
             
            
        m_=LL_-2    
        c2_=(l_+m_+2d0)*(l_-m_-1d0)&
              /dble(l_-m_)/(l_+m_+1d0)
            coef2_=dsqrt(c2_)
         c(l_+1)=  coef2_ 
        
        f=0d0
        f(0+1)=-pl0_(l_+1)
        f(l_-1)=-p_(l_+1)*c(l_+1)
        c(1)=c(1)/b(1)
        f(1)=f(1)/b(1)
        do  m_=1,l_-2
            c1_=4d0*(dble(m_)+1.d0)**2&
                 /dble(l_-m_)/(l_+m_+1d0)
            coef1_=dsqrt(c1_)

            c2_=(l_+m_+2d0)*(l_-m_-1d0)&
              /dble(l_-m_)/(l_+m_+1d0)
            coef2_=dsqrt(c2_)
            
            b(m_+1)=coef1_*x_/y_
            c(m_+1)=coef2_
            
  
         c(m_+1)=c(m_+1)/(b(m_+1)-c(m_))
         
         f(m_+1)=(f(m_+1)-f(m_))/(b(m_+1)-c(m_))
         
         
   
        end do
        p_(l_)=f(l_-1)

        p_(1)=pl0_(l_+1)
        
        
        do m_=l_-3,0,-1
         p_(m_+2)=f(m_+1)-c(m_+1)*p_(m_+3)
        end do
        
    else 
    
    


    do m_=1,LL_
        coef_=-dsqrt((2.d0*dble(m_)+1.d0)/(2.d0*dble(m_)))

        p_(m_+1)=coef_*p_(m_-1+1)*y_
            !print*,x(j)

    end do
   
    m_=LL_
        coef_=sqrt(2d0*m_)
        !print*,coef_

        p_(m_-1+1)=-coef_*x_/y_*p_(m_+1)

    

    l_=LL_
    
    
        do  m_=l_-2,0,-1

            c1_=4d0*(dble(m_)+1.d0)**2&
                 /dble(l_-m_)/(l_+m_+1d0)
            coef1_=dsqrt(c1_)

            c2_=(l_+m_+2d0)*(l_-m_-1d0)&
              /dble(l_-m_)/(l_+m_+1d0)
            coef2_=dsqrt(c2_)


            p_(m_+1)=-coef1_*x_/y_*p_(m_+1+1)-coef2_*p_(m_+2+1)
            
        end do
    
    end if
    
    
end subroutine


subroutine findsingu(i,j,MM,q00,q10,q01,q11,u00,u10,u01,u11,&
 qx00,qx10,qx01,qx11,ux00,ux10,ux01,ux11,&
 qy00,qy10,qy01,qy11,uy00,uy10,uy01,uy11,&
 qxnp,qynp,uxnp,uynp,tnpf,pnpf,nflag)
     implicit none

         integer, intent(in) :: i, j, MM
         integer, intent(out) :: nflag
         doubleprecision, intent(in) :: q00,q10,q01,q11,u00,u10,u01,u11,&
 qx00,qx10,qx01,qx11,ux00,ux10,ux01,ux11,&
 qy00,qy10,qy01,qy11,uy00,uy10,uy01,uy11
         doubleprecision, intent(out) :: qxnp,qynp,uxnp,uynp,tnpf,pnpf

         doubleprecision xq(5),yq(5),zq(5),xu(5),yu(5),zu(5)
         doubleprecision qx(5),qy(5),ux(5),uy(5)
         doubleprecision a, a1, a2, ab, ac, alpha, ax, ay, az
         doubleprecision bc, betta, bx, by, bz, c
         doubleprecision ca, ca1, ca2, clow, ctnp
         doubleprecision cup, cx, cy, cz, delta, den
         doubleprecision  f, f1, f2, h
         doubleprecision  p0, p1, pl, plow, pnp, pr
         integer nq, nu, kl, kll
         doubleprecision psi, pup, qxnp1, qxnp2, qynp1, qynp2
         doubleprecision r, s, t0, t1, tl, tlow, tlower, tnp
         doubleprecision tr, tup, tupper
         doubleprecision uxnp1, uxnp2, uynp1, uynp2
         doubleprecision x, x1, x2, xnp
         doubleprecision y, y1, y2, ynp
         doubleprecision z, z1, z2, znp
!c         dimension F(0:4098,0:1026)
!c         dimension xf1(0:10000),yf1(0:10000)
         
         
            pi=4.d0*datan(1.d0)
            
!c---------------------input----------------------

            
 !         print 1000
 !1000     format('input i (phi), j (Theta):',',' )
 !         read *, i,j


!c-------------number of pixels along theta--------
!c--------------1<=j<=MM-2-------------------------
!c--------------0<=i<=2*MM-1-----------------------          
!            MM=180
            h=pi/dble(MM)
            p0=dble(i)*h
            t0=dble(j)*h
            p1=p0+h
            t1=t0+h


            
!c-----------------------------------------------            
          nq=0
          nu=0
          nflag=0

!c-----------------along phi---------------------
!c----------------- q lower--------------------------
          
          if (q00.le.0.d0.and.q10.gt.0.d0.or.&
             q00.ge.0.d0.and.q10.lt.0.d0) then
             nq=nq+1

             x1=dsin(t0)*dcos(p0)
             y1=dsin(t0)*dsin(p0)
             z1=dcos(t0)
             x2=dsin(t0)*dcos(p1)
             y2=dsin(t0)*dsin(p1)
             z2=dcos(t0)
             ab=x1*x2+y1*y2+z1*z2
             psi=dacos(ab)
             delta=psi*q00/(q00-q10)
             ac=dcos(delta)
             bc=dcos(psi-delta)
             den=1.d0-ab*ab
             alpha=(ac-ab*bc)/den
             betta=(bc-ab*ac)/den
             x=alpha*x1+betta*x2
             y=alpha*y1+betta*y2
             z=alpha*z1+betta*z2
             r=dsqrt(x*x+y*y+z*z)
             
             xq(nq)=x/r
             yq(nq)=y/r
             zq(nq)=z/r

!c-----------------interpolation-----------------------
  
             f1=qx00
             f2=qx10
             f=(f2-f1)/psi*delta+f1
             qx(nq)=f
             
             f1=qy00
             f2=qy10
             f=(f2-f1)/psi*delta+f1
             qy(nq)=f
             
!c-----------------------------------------------------             
            else
            endif

!c------------------------q upper---------------------            
            
           if (q01.le.0.d0.and.q11.gt.0.d0.or.&
              q01.ge.0.d0.and.q11.lt.0.d0) then
              nq=nq+1

             x1=dsin(t1)*dcos(p0)
             y1=dsin(t1)*dsin(p0)
             z1=dcos(t1)
             x2=dsin(t1)*dcos(p1)
             y2=dsin(t1)*dsin(p1)
             z2=dcos(t1)
             ab=x1*x2+y1*y2+z1*z2
             psi=dacos(ab)
             delta=psi*q01/(q01-q11)
             ac=dcos(delta)
             bc=dcos(psi-delta)
             den=1.d0-ab*ab
             alpha=(ac-ab*bc)/den
             betta=(bc-ab*ac)/den
             x=alpha*x1+betta*x2
             y=alpha*y1+betta*y2
             z=alpha*z1+betta*z2
             r=dsqrt(x*x+y*y+z*z)
             
             xq(nq)=x/r
             yq(nq)=y/r
             zq(nq)=z/r


!c-----------------interpolation-----------------------
  
             f1=qx01
             f2=qx11
             f=(f2-f1)/psi*delta+f1
             qx(nq)=f
             
             f1=qy01
             f2=qy11
             f=(f2-f1)/psi*delta+f1
             qy(nq)=f

!c--------------------------------------------------------
             
            else
            endif

!c---------------------u lower------------------------------
            
            if (u00.le.0.d0.and.u10.gt.0.d0.or.&
               u00.ge.0.d0.and.u10.lt.0.d0) then
             nu=nu+1

             x1=dsin(t0)*dcos(p0)
             y1=dsin(t0)*dsin(p0)
             z1=dcos(t0)
             x2=dsin(t0)*dcos(p1)
             y2=dsin(t0)*dsin(p1)
             z2=dcos(t0)
             ab=x1*x2+y1*y2+z1*z2
             psi=dacos(ab)
             delta=psi*u00/(u00-u10)
             ac=dcos(delta)
             bc=dcos(psi-delta)
             den=1.d0-ab*ab
             alpha=(ac-ab*bc)/den
             betta=(bc-ab*ac)/den
             x=alpha*x1+betta*x2
             y=alpha*y1+betta*y2
             z=alpha*z1+betta*z2
             r=dsqrt(x*x+y*y+z*z)
             
             xu(nu)=x/r
             yu(nu)=y/r
             zu(nu)=z/r

!c-----------------interpolation-----------------------
  
             f1=ux00
             f2=ux10
             f=(f2-f1)/psi*delta+f1
             ux(nu)=f
             
             f1=uy00
             f2=uy10
             f=(f2-f1)/psi*delta+f1
             uy(nu)=f

!c--------------------------------------------------------
             
            else
            endif

!c-------------------u upper----------------            
            
           if (u01.le.0.d0.and.u11.gt.0.d0.or.&
              u01.ge.0.d0.and.u11.lt.0.d0) then
              nu=nu+1

             x1=dsin(t1)*dcos(p0)
             y1=dsin(t1)*dsin(p0)
             z1=dcos(t1)
             x2=dsin(t1)*dcos(p1)
             y2=dsin(t1)*dsin(p1)
             z2=dcos(t1)
             ab=x1*x2+y1*y2+z1*z2
             psi=dacos(ab)
             delta=psi*u01/(u01-u11)
             ac=dcos(delta)
             bc=dcos(psi-delta)
             den=1.d0-ab*ab
             alpha=(ac-ab*bc)/den
             betta=(bc-ab*ac)/den
             x=alpha*x1+betta*x2
             y=alpha*y1+betta*y2
             z=alpha*z1+betta*z2
             r=dsqrt(x*x+y*y+z*z)
             
             xu(nu)=x/r
             yu(nu)=y/r
             zu(nu)=z/r

!c-----------------interpolation-----------------------
  
             f1=ux01
             f2=ux11
             f=(f2-f1)/psi*delta+f1
             ux(nu)=f
             
             f1=uy01
             f2=uy11
             f=(f2-f1)/psi*delta+f1
             uy(nu)=f

!c--------------------------------------------------------

             
            else
            endif            
            
!c-------------------end of along phi---------------


!c------------------along theta----------------------

!c-----------------q left---------------------------
            
             if (q00.le.0.d0.and.q01.gt.0.d0.or.&
                q00.ge.0.d0.and.q01.lt.0.d0) then
                nq=nq+1
                
                delta=h*q00/(q00-q01)
                tl=t0+delta
                pl=p0

                xq(nq)=dsin(tl)*dcos(pl)
                yq(nq)=dsin(tl)*dsin(pl)
                zq(nq)=dcos(tl)

!c-----------------interpolation-----------------------
  
             f1=qx00
             f2=qx01
             f=(f2-f1)/h*delta+f1
             qx(nq)=f
             
             f1=qy00
             f2=qy01
             f=(f2-f1)/h*delta+f1
             qy(nq)=f

             
!c--------------------------------------------------------
                
                  else
              endif

!c-----------------q right---------------------------
            
             if (q10.le.0.d0.and.q11.gt.0.d0.or.&
                q10.ge.0.d0.and.q11.lt.0.d0) then
                nq=nq+1
                
                delta=h*q10/(q10-q11)
                tr=t0+delta
                pr=p1

                xq(nq)=dsin(tr)*dcos(pr)
                yq(nq)=dsin(tr)*dsin(pr)
                zq(nq)=dcos(tr)

!c-----------------interpolation-----------------------
  
             f1=qx10
             f2=qx11
             f=(f2-f1)/h*delta+f1
             qx(nq)=f
             
             f1=qy10
             f2=qy11
             f=(f2-f1)/h*delta+f1
             qy(nq)=f

!c--------------------------------------------------------
                
             else
             endif

!c-----------------u left---------------------------
            
             if (u00.le.0.d0.and.u01.gt.0.d0.or.&
                u00.ge.0.d0.and.u01.lt.0.d0) then
                nu=nu+1
                
                delta=h*u00/(u00-u01)
                tl=t0+delta
                pl=p0

                xu(nu)=dsin(tl)*dcos(pl)
                yu(nu)=dsin(tl)*dsin(pl)
                zu(nu)=dcos(tl)

!c-----------------interpolation-----------------------
  
             f1=ux00
             f2=ux01
             f=(f2-f1)/h*delta+f1
             ux(nu)=f
             
             f1=uy00
             f2=uy01
             f=(f2-f1)/h*delta+f1
             uy(nu)=f

!c--------------------------------------------------------
                
                  else
              endif

!c-----------------u right---------------------------
            
             if (u10.le.0.d0.and.u11.gt.0.d0.or.&
                u10.ge.0.d0.and.u11.lt.0.d0) then
                nu=nu+1
                
                delta=h*u10/(u10-u11)
                tr=t0+delta
                pr=p1

                xu(nu)=dsin(tr)*dcos(pr)
                yu(nu)=dsin(tr)*dsin(pr)
                zu(nu)=dcos(tr)

!c-----------------interpolation-----------------------
  
             f1=ux10
             f2=ux11
             f=(f2-f1)/h*delta+f1
             ux(nu)=f
             
             f1=uy10
             f2=uy11
             f=(f2-f1)/h*delta+f1
             uy(nu)=f

!c--------------------------------------------------------
                
             else
                endif

             
             

            if (nq.eq.2.and.nu.eq.2) then

!c------------------------------------------------------
!c------------------intersection------------------------
               
               x1=xq(1)
               y1=yq(1)
               z1=zq(1)

               x2=xq(2)
               y2=yq(2)
               z2=zq(2)
               
               ax=y1*z2-z1*y2
               ay=-x1*z2+z1*x2
               az=x1*y2-y1*x2

               x1=xu(1)
               y1=yu(1)
               z1=zu(1)

               x2=xu(2)
               y2=yu(2)
               z2=zu(2)
               
               bx=y1*z2-z1*y2
               by=-x1*z2+z1*x2
               bz=x1*y2-y1*x2

               cx=ay*bz-az*by
               cy=-ax*bz+az*bx
               cz=ax*by-ay*bx
               r=dsqrt(cx*cx+cy*cy+cz*cz)
!c--------------------------------------------               
               cx=cx/r
               cy=cy/r
               cz=cz/r
               
               do 1 kl=1,2

                  cx=-cx
                  cy=-cy
                  cz=-cz
                                 
               ctnp=cz
               tnp=dacos(ctnp)
                              
               s=cy/dsin(tnp)
               c=cx/dsin(tnp)

               if (s.ge.1.d0) then
                  s=1.d0
               else
               endif

               if (s.le.-1.d0) then
                  s=-1.d0
               else
               endif

               if (c.ge.1.d0) then
                  c=1.d0
               else
               endif

               if (c.le.-1.d0) then
                  c=-1.d0
               else
                endif 
              
               if (s.ge.0) then                  
                 pnp=dacos(c)
                   else
                 pnp=2.d0*pi-dacos(c)
               endif               
!c------------------------------------------------------               

                if (pnp.ge.p0.and.pnp.lt.p1) then                   

                  tnpf=tnp
                  pnpf=pnp
                  
                  xnp=dsin(tnpf)*dcos(pnpf)
                  ynp=dsin(tnpf)*dsin(pnpf)
                  znp=dcos(tnpf)
                  
!c-----------------------checking lower limit------------
                  
             x1=dsin(t0)*dcos(p0)
             y1=dsin(t0)*dsin(p0)
             z1=dcos(t0)
             x2=dsin(t0)*dcos(p1)
             y2=dsin(t0)*dsin(p1)
             z2=dcos(t0)
             ax=y1*z2-z1*y2
             ay=-x1*z2+z1*x2
             az=x1*y2-y1*x2

             x1=dsin(tnp)*dcos(pnp)
             y1=dsin(tnp)*dsin(pnp)
             z1=dcos(tnp)
             x2=0.d0
             y2=0.d0
             z2=1.d0
             bx=y1*z2-z1*y2
             by=-x1*z2+z1*x2
             bz=x1*y2-y1*x2
             cx=ay*bz-az*by
             cy=-ax*bz+az*bx
             cz=ax*by-ay*bx
             r=dsqrt(cx*cx+cy*cy+cz*cz)

             cx=cx/r
             cy=cy/r
             cz=cz/r

              do 2 kll=1,2

                  cx=-cx
                  cy=-cy
                  cz=-cz
                                 
               clow=cz
               tlow=dacos(clow)
                              
               s=cy/dsin(tlow)
               c=cx/dsin(tlow)

               if (s.ge.1.d0) then
                  s=1.d0
               else
               endif

               if (s.le.-1.d0) then
                  s=-1.d0
               else
               endif

               if (c.ge.1.d0) then
                  c=1.d0
               else
               endif

               if (c.le.-1.d0) then
                  c=-1.d0
               else
                endif 
              
               if (s.ge.0) then                  
                 plow=dacos(c)
                   else
                 plow=2.d0*pi-dacos(c)
               endif        

               if (plow.ge.p0.and.plow.lt.p1) then
!c              print *, kl,kll,tlow/pi*180.d0,plow/pi*180.d0
               tlower=tlow   
               else
                endif  
                  

 2           continue
             
!c-----------------------checking upper limit------------
             
             x1=dsin(t1)*dcos(p0)
             y1=dsin(t1)*dsin(p0)
             z1=dcos(t1)
             x2=dsin(t1)*dcos(p1)
             y2=dsin(t1)*dsin(p1)
             z2=dcos(t1)
             ax=y1*z2-z1*y2
             ay=-x1*z2+z1*x2
             az=x1*y2-y1*x2

             x1=dsin(tnp)*dcos(pnp)
             y1=dsin(tnp)*dsin(pnp)
             z1=dcos(tnp)
             x2=0.d0
             y2=0.d0
             z2=1.d0
             bx=y1*z2-z1*y2
             by=-x1*z2+z1*x2
             bz=x1*y2-y1*x2
             cx=ay*bz-az*by
             cy=-ax*bz+az*bx
             cz=ax*by-ay*bx
             r=dsqrt(cx*cx+cy*cy+cz*cz)

             cx=cx/r
             cy=cy/r
             cz=cz/r

              do 3 kll=1,2

                  cx=-cx
                  cy=-cy
                  cz=-cz
                                 
               cup=cz
               tup=dacos(cup)
                              
               s=cy/dsin(tup)
               c=cx/dsin(tup)

               if (s.ge.1.d0) then
                  s=1.d0
               else
               endif

               if (s.le.-1.d0) then
                  s=-1.d0
               else
               endif

               if (c.ge.1.d0) then
                  c=1.d0
               else
               endif

               if (c.le.-1.d0) then
                  c=-1.d0
               else
                endif 
              
               if (s.ge.0) then                  
                 pup=dacos(c)
                   else
                 pup=2.d0*pi-dacos(c)
               endif        

               if (pup.ge.p0.and.pup.lt.p1) then
                tupper=tup  
               else
                endif  
                  

 3           continue
             


               if (tnpf.ge.tlower.and.tnpf.le.tupper) then

!c-----------------q interpolation--------------------

               ca1=xnp*xq(1)+ynp*yq(1)+znp*zq(1)
                if (ca1.ge.1.d0) then
                   ca1=1.d0
                else
                endif
                a1=dacos(ca1)
                
               ca2=xnp*xq(2)+ynp*yq(2)+znp*zq(2)
                 if (ca2.ge.1.d0) then
                   ca2=1.d0
                else
                   endif
                   a2=dacos(ca2)

                  ca=xq(1)*xq(2)+yq(1)*yq(2)+zq(1)*zq(2)
                 if (ca.ge.1.d0) then
                   ca=1.d0
                else
                   endif
                   a=dacos(ca)  

                   f1=qx(1)
                   f2=qx(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   qxnp1=f

                   f1=qy(1)
                   f2=qy(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   qynp1=f

                   f1=ux(1)
                   f2=ux(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   uxnp1=f

                   f1=uy(1)
                   f2=uy(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   uynp1=f


!c-----------------u interpolation--------------------

               ca1=xnp*xu(1)+ynp*yu(1)+znp*zu(1)
                if (ca1.ge.1.d0) then
                   ca1=1.d0
                else
                endif
                a1=dacos(ca1)
                
               ca2=xnp*xu(2)+ynp*yu(2)+znp*zu(2)
                 if (ca2.ge.1.d0) then
                   ca2=1.d0
                else
                   endif
                   a2=dacos(ca2)

                  ca=xu(1)*xu(2)+yu(1)*yu(2)+zu(1)*zu(2)
                 if (ca.ge.1.d0) then
                   ca=1.d0
                else
                   endif
                   a=dacos(ca)  

                   f1=qx(1)
                   f2=qx(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   qxnp2=f

                   f1=qy(1)
                   f2=qy(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   qynp2=f

                   f1=ux(1)
                   f2=ux(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   uxnp2=f

                   f1=uy(1)
                   f2=uy(2)
                   f=(f2-f1)/(a1+a2)*a1+f1
                   uynp2=f  
                   
                   qxnp=(qxnp1+qxnp2)/2.d0
                   qynp=(qynp1+qynp2)/2.d0
                   uxnp=(uxnp1+uxnp2)/2.d0
                   uynp=(uynp1+uynp2)/2.d0                   

                   nflag=1
                   
!c-----------------------------------------------------------                   

                    else
                 endif

                 else
              endif
                  
 1             continue
               
                    else
                       
               endif

!c---------------------output-------------------------------

              if (nflag.eq.1) then    
               
               tnpf=tnpf/h
               pnpf=pnpf/h
            else
            !   print *, nflag
             endif  
!c-----------------------------------------------------------
               
            

end subroutine             





end program polar
