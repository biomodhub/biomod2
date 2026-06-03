C Output from Public domain Ratfor, version 1.03
      subroutine fexit(msg)
      character*(*) msg
      nc = len(msg)
      call fexitc(msg, nc)
      end
