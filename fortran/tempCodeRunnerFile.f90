program AtCoder
    implicit none
    character(len = 10**5) :: s
    integer,dimension(0:10**5) :: flag = 0
    integer :: i
    read *,s
    flag(0) = 1
    do i = 0,len_trim(s)
        if(flag(i) == 1)then
            if(s(i+1:i+5) == "dream")then
                flag(i+5) = 1
            end if
            if(s(i+1:i+7) == "dreamer")then
                flag(i+7) = 1
            end if
            if(s(i+1:i+5) == "erase")then
                flag(i+5) = 1
            end if
            if(s(i+1:i+6) == "eraser")then
                flag(i+6) = 1
            end if
        end if
    end do
    if(flag(len_trim(s)) == 1)then
        print *,"YES"
    else
        print *,"NO"
    end if
contains
    
end program AtCoder