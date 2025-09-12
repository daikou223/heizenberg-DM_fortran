program mask
    integer::n,a,b,minT,maxT,maxMask,minMask,restMask
    n = 30
    a = 0
    b = 9
    maxT = max(a,b)
    minT = min(a,b)
    restMask = (ishft(1,minT))-1
    minMask = (ishft(1,maxT))-1-restMask
    maxMask = (ishft(1,N))-1-minMask-restmask
    print *,restMask
    print *,minMask
    print *,maxMask
end program mask