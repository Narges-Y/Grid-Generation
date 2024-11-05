program Grid_generation
    implicit none

    real,parameter::pi= 3.14159;
    real :: x(61,31) , y(61,31) , xc(31) , yc(31)
    real :: r, dt, dx, dy
    integer :: pn, pm, i, j, imax, jmax

    open(19,file='grid.plt')

    !chord discretization

    pn = 30
    pm = 30
    imax = 2*pn+1
    jmax = pm+1
    r = 4.0
    dt = pi/pn

    !boundary discretization

    call meth(0.005 , 0.02 , 1.0 , pn+1 , xc)
    call airfoil(0.12 , 1.0 , pn+1 , xc , yc)


    do i=1,pn+1
        x(i,1) = xc(i)-0.5
        y(i,1) = yc(i)

        end do

        do i=pn+2,imax
            x(i,1)=xc(2*pn-i+2)-0.5
            y(i,1)=-yc(2*pn-i+2)
            end do

            !O-type

            do i=1,pn+1
            x(i,jmax)=-r*cos((i-1)*dt)
            y(i,jmax)=r*sin((i-1)*dt)
        end do

        do i=pn+2,imax
            x(i,jmax)=x(2*pn+2-i,jmax)
            y(i,jmax)=-y(2*pn+2-i,jmax)
        end do




        do i=1,imax
            dx=(((x(i,jmax))-x(i,1))/pm)
            dy=(((y(i,jmax))-y(i,1))/pm)

            do j=2,jmax-1
                x(i,j)=x(i,1)+(j-1)*dx
                y(i,j)=y(i,1)+(j-1)*dy
            end do
        end do

        write(19,*)'variables = x , y ,z'
        write(19,*)'zone i= ' , imax , 'j= ' ,jmax

        do j=1,jmax
        do i=1,imax
            write(19,*) x(i,j) , y(i,j) , '0.0'
        end do
        end do


            !Method
            contains

            subroutine meth(dsa,dsb,smax,jmax,u)
              implicit none

            real:: dsa , dsb , smax , am , bm , Tt , a , b , c , etha , u(jmax) , delta , al , bl
            integer :: jmax , i , j
            am = sqrt(dsb/dsa)
            bm = smax/((jmax-1)*sqrt(dsa*dsb))
            Tt=0.0001
            a=0.01
            b=10.0


                do i=1,100
                    c=(a+b)/2
                    al=(sinh(c)/c)-bm
                    bl=(sinh(a)/a)-bm

                    if(al==0 .or. (b-a)/2<Tt) exit
                    if(al*bl >= 0) then
                        a=c
                    else
                        b=c
                    end if
                    end do

            delta=c
            do j=1,jmax
                etha=real(j-1)/((jmax)-1)
                u(j)=0.5+tanh(delta*(etha-0.5))/(2.0*tanh(delta/2.0))
                end do

            end subroutine


        subroutine airfoil(t , c , n , x , y)
            integer :: i , n
            real :: t , c , y(n) , x(n)

            do i=1,n
                y(i)= 5.0*t*c*(0.2969*sqrt (x(i)/c)+ (((-0.1015)*(x(i)/c)+0.2843)*(x(i)/c)-0.3516)*(x(i)/c)-0.1260*(x(i)/c))

                if(x(i)==1)then
                    y(i)=0
                else if(x(i)==0)then
                    y(i)=0
                end if
                end do
            end subroutine

end program

