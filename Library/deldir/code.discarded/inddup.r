subroutine inddup(x,y,n,rw,frac,dup)
implicit double precision(a-h,o-z)
logical dup(n)
dimension x(n), y(n), rw(4)

xtol = frac*(rw(2)-rw(1))
ytol = frac*(rw(4)-rw(3))

dup(1) = .false.
do i = 2,n {
	dup(i) = .false.
	do j = 1,i-1 {
		dx = abs(x(i)-x(j))
		dy = abs(y(i)-y(j))
		if(dx < xtol & dy < ytol) {
			dup(i) = .true.
			break
		}
	}
}

return
end
