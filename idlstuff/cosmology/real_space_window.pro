;; F.T. of real-space window function for calculating sigma(M)
function real_space_window, x
	W = (3./(x^3)) * (sin(x) - x*cos(x))
	return, W
end
