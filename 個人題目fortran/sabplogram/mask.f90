program mask
    integer::N,mask,state,staetInd,
    mask = (ishft(1,N))-1
    if(idToInd(ishft(state,-2)) <= mask) then
        stateInd = idToInd(ishft(state,-2))
    else
        if(state == states(iand(idToInd(ishft(state,-2)),mask))) then
            stateInd = iand(idToInd(ishft(state,-2)),mask)
        end if
        if(state == states(ishft(idToInd(ishft(state,-2)),-N)))then
            stateInd =  ishft(idToInd(ishft(state,-2)),-N)
        end if
    end if
    !fripstateのほう
    if(idToInd(ishft(fripState,-2)) <= mask)then
        fripStateInd = idToInd(ishft(fripState,-2))
    else
        if(fripState == states(iand(idToInd(ishft(fripState,-2)),mask))) then
            fripStateInd = iand(idToInd(ishft(fripState,-2)),mask)
        end if
        if(fripState == states(ishft(idToInd(ishft(fripState,-2)),-N))) then
            fripStateInd =  ishft(idToInd(ishft(fripState,-2)),-N)
        end if
    end if
    if(stateInd < 1 .or. stateInd > ALL_STATE_NUM .or. fripStateInd < 1 .or. &
        fripStateInd > ALL_STATE_NUM) then
        print *,stateInd,fripStateInd,state,fripState
    end if
end program mask