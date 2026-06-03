subroutine fexit(msg)
character*(*) msg
nc = len(msg)
call fexitc(msg, nc)
end
