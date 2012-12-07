module stimeparam 


type :: time

logical :: acceleration_scheme, velocity_scheme
integer :: ntimeMax, NtimeMin, nSnap, ntrace, ncheck
real :: alpha, beta, gamma, duration, Time_snapshots, dtmin, rtime

end type

end module stimeparam
