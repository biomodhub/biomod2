
test_sequential <- function(max=100, nb=1000, display_progress=TRUE) {
	.Call("test_sequential_wrapper2", max, nb, display_progress, PACKAGE = "RcppProgressETA")
	invisible()
}

test_multithreaded <- function(max=100, nb=1000, threads=0, display_progress=TRUE) {
	.Call("test_multithreaded_wrapper2", max, nb, threads, display_progress, PACKAGE = "RcppProgressETA")
	invisible()
}
