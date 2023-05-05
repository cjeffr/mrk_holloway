! Program MRK HOLLOWAY 10/29/91
! Program Unit MRKHOL Entry:00000444 Options: ILSDNB

program mrk_holloway
    ! Program takes input in the form of user input and calculates an
    ! Isochore from temperature (deg c) and molar volume of constituents
    implicit none

    ! Program Arguments
    integer, parameter :: dp = selected_real_kind(15, 307)
    character(len=80) :: zout
    character(len=5), dimension(8):: names=(/' CO2', '  CO', ' CH4', '  H2', &
                                          ' H20', ' H2S', ' SO2', '  N2'/)
    real(kind=dp), dimension(8) :: x

    ! Local Storage
    real(kind=dp), dimension(8) :: formwt=(/44.d0, 28.d0, 16.d0, 2.d0, &
                                            18.d0, 34.d0, 64.06d0, 28.0d0/)
    real(kind=dp), dimension(5) :: voli, pout, den
    real(kind=dp) :: tstart, vstart, form, tk, pstart, t, tc, asum, bsum
    real(kind=dp) :: vol, p, aterm, bterm
    integer :: mixnum = 8, i, inext, j
    logical :: is_one

    call run_all()

        
    contains
    subroutine additional_runs()
        print *, 'Do you want to do another composition? No=0, Yes=1'
        read(*,*) inext
        if (inext.eq.0) then
            stop
        end if
        if (inext.eq.1) then
            call run_all()
            
        end if
    end subroutine additional_runs

    subroutine run_all()
        call obtain_input(x, names, zout, tstart, vstart)
        ! Check for xtotal = 1
    call check_totals(x, is_one)
    if (is_one) then
        continue
    else
        print *, 'The mole fractions do not add up to 1.0! Start again'
        call obtain_input(x, names, zout, tstart, vstart)
    end if

    print *,' The mole fractions are:'
    print '(8(5X,A4))', names
    print '(8(f9.3))', x
    print '(a,f6.0,a)', 'The starting temperatature = ', tstart, ' deg. c'
    print '(a,f6.3,a)', 'The starting molar volume = ', vstart, ' cm^3/mol'


    call calc_form_weight(form, x, mixnum)

   
    
    tk = tstart + 273.15
    ! pstart = 1.d3
    call calc_isochore_volume(form, vstart, den, voli)
    print '(30x, a)', 'Isochore pressures in bars'
    print '(a, 5f10.3)', &
        ' T deg. c mol vol =',voli 
    print '(10x,a, 5f10.4)', 'density =', den


    call temp_pressure_calc(x, tk, t, pout, tc, voli)
    call additional_runs()

end subroutine run_all
    subroutine check_totals(x, is_one)
        implicit none
        ! Input
        real(kind=dp), dimension(8), intent(in) :: x
        ! Local
        real(kind=dp) :: xtotal
        logical, intent(out) :: is_one
        xtotal = sum(x)
        if (xtotal.eq.1.00) then
            is_one = .True.
        else
            is_one = .False.
        end if
    end subroutine check_totals

    subroutine temp_pressure_calc(x, tk, t, pout, tc, voli)
        implicit none
        ! Imports
        real(kind=dp), dimension(8), intent(in) :: x
        real(kind=dp), intent(in) :: tk
        real(kind=dp), intent(in), dimension(5) :: voli
        
       

        ! Exports
        real(kind=dp), intent(out), dimension(5) :: pout
        real(kind=dp), intent(out) :: tc
        real(kind=dp), intent(out) :: t

        ! Local
        integer :: i, j
        real(kind=dp) :: vol, p
        real(kind=dp) :: rbar = 83.117d0
        real(kind=dp) :: asum, bsum, aterm, bterm

        do i=1,11
            t = tk - 1.d2 + 1.d2 * i
            tc = t - 273.15
            call mrkmix(tk, x, bsum, asum)
            ! Calculate isochore
            do j=1,5
                vol = voli(j)
                p = 0.d0
                if (vol.lt.(bsum - 1.0d0) ) then
                    continue
                else 
                   
                    aterm = ((rbar * t) / (vol - bsum))
                    bterm = asum / ((dsqrt(t)) * ((vol * vol) + (bsum * vol)))
                    p = aterm - bterm
                pout(j) = p
                    
                end if
            end do
            print '(f10.0, 10x, 5f10.0)', tc, pout
        end do

    end subroutine temp_pressure_calc

    subroutine calc_isochore_volume(form, vstart, den, voli)
        implicit none
        ! Import
        real(kind=dp), intent(in) :: form
        real(kind=dp), intent(in) :: vstart
        ! Export
        real(kind=dp), dimension(5), intent(out) :: den, voli

        ! Calculate total volume for an isochore
        do i=1,5
            voli(i) = vstart + 1.d1 * i - 1.d1
            den(i) = form/voli(i)
        end do
    end subroutine calc_isochore_volume

    subroutine obtain_input(x, names, zout, tstart, vstart)
        implicit none
        character(len=80), intent(out) :: zout
        character, dimension(8), intent(in) :: names
        real(kind=dp), dimension(8), intent(out) :: x
        real(kind=dp), intent(out) :: tstart, vstart

        print *, 'Enter name of output file'
        read(*, *) zout
        open(unit=8, file=zout)
        print *, 'Enter T (degC) and Molar Volume (cc/mole)'
        read(*,*) tstart, vstart
        if(tstart .LT. 1.d-2) then
            stop
        end if
        print *, 'Enter mole fraction CO2'
        read(*,*) x(1)
        print *, 'Enter mole fraction CO'
        read(*,*) x(2)
        print *, 'Enter mole fraction CH4'
        read(*,*) x(3)
        print *, 'Enter mole fraction H2'
        read(*,*) x(4)
        print *, 'Enter mole fraction H2O'
        read(*,*) x(5)
        print *, 'Enter mole fraction H2S'
        read(*,*) x(6)
        print *, 'Enter mole fraction SO2'
        read(*,*) x(7)
        print *, 'Enter mole fraction N2'
        read(*,*) x(8)

    end subroutine obtain_input


    subroutine calc_form_weight(form, x, mixnum)
        implicit none
        ! Import
        integer, intent(in) :: mixnum
        real(kind=dp), dimension(8), intent(in) :: x
        ! export
        real(kind=dp), intent(out) :: form

        ! Calculate mean formula weight
        form = 0.0
        do i=1,mixnum
            form = form + x(i) * formwt(i)
        end do

    end subroutine calc_form_weight

    subroutine mrkmix(t_kelvin, y, bsum, asum)
        implicit none
        ! Import
        real(kind=dp), dimension(8), intent(in) :: y
        real(kind=dp), intent(in) :: t_kelvin

        ! Export
        real(kind=dp), intent(out) :: bsum, asum

        ! Locals
        real, dimension(8) :: a = (/46.d6, 16.98d6, 31.59d6, 3.56d6, &
                                    35.d6, 87.9d6, 142.6d6, 15.382d6/)
        real, dimension(8) :: b = (/2.97d1, 2.738d1, 2.9703d1, 1.515d1, &
                                   1.46d1, 2.d1, 3.94d1, 2.68d1/)
        real(kind=dp) :: tcelsius, r2t, rt, ah20m
        real(kind=dp) :: aco2m, co2h20, xk, t, r=82.05d0
        integer :: mixnum=8, i, j

        if (t_kelvin.lt.1.d-4) then
            t = 1.d0
        else
            t = t_kelvin
        end if

        tcelsius = t - 273.15
        r2t = r * r * t**2.5
        rt = r * t**1.5
        ah20m = 166.8 - .19308 * tcelsius + .1864d-3 * tcelsius * tcelsius - .71288d-7 * tcelsius**3
        if (tcelsius.lt.6.d2) then
            ah20m = 4.221d3 - 3.1227d1 * tcelsius + 8.7485d-2 *  tcelsius**2 - &
            1.07295d-4 * tcelsius**3 + 4.86111d-8 * tcelsius**4
            

        else if (tcelsius.gt.1200) then
            ah20m = 140. - 0.050 * tcelsius
        end if
        
        ah20m = ah20m * 10.d5
        aco2m = 73.03 - 0.0714 * tcelsius + 2.157d-5 * tcelsius**2
        aco2m = aco2m * 10.d5
        xk = exp(-11.071 + ( 5953./t) - (2.746d6/(t**2)) + &
            (4.646d8/(t**3)))
        co2h20 = xk * 0.5 * r2t
        co2h20 = co2h20 + sqrt(a(1) * a(5))
        
        asum = 0.0
        bsum = 0.0

        do i = 1,mixnum
            bsum = bsum + b(i) * y(i)
            do j=1,mixnum
                
                if (i.eq.j) then
                    if ((i.ne.5) .and. (i.ne.1)) then
                        asum = asum + y(i) * y(j) * a(i)

                    else if (i.eq.5) then
                        asum = asum + y(i) * y(j) * ah20m
                        
                    else if (i.eq.1) then
                        asum = asum + y(i) * y(j) * aco2m  
                    end if
                

                else if (((i.eq.5) .and. (j.eq.1)) .or. ((i.eq.1) .and.(j.eq.5))) then
                    asum = asum + y(i) * y(j) * co2h20
                else 
                    asum = asum + y(i) * y(j) * sqrt(a(i)*a(j))
                    
                end if 
        
            end do
        end do
        asum = asum/1.013
            
    end subroutine mrkmix
end program mrk_holloway