using DSP

function conv_py(times, noise)
	py"""
	import numpy as np
	import sys
	sys.path.append("src/")
	import integration as itg
	def conv__(t, noise):
		t__ = np.array(t)
		noise__ = np.array(noise)
		def linear(t):
			return t
		conv = itg.convolution(linear, noise__, t__)
		return conv
	"""
    cnv = py"conv__"(times, noise)
    return cnv
end

noise = get_noise(h, noise_steps*N, final_time, 1, DATA_PATH)
times = [i*final_time/N for i in 0:N-1]

cnv1 = conv_py(times, noise)
cnv2 = conv(times, noise)
plot(cnv1)
plot!(cnv2)


using MittagLeffler

h = 0.1
fd_order = 2 - 2*h
N = 100
final_time = 3

times = [i*final_time/N for i in 0:N-1]
anl = []
for t in times
	println(t)
	append!(anl, [mittleff(2, 1, t)])
end

function mtlf_py(times, noise)
	py"""
	import numpy as np
	import sys
	sys.path.append("src/")
	import mittag_leffler as ml
	def mtlf(t):
		t__ = np.array(t)
		return ml.mittag_leffler_vector(t__, 2-order, 2)
	"""
    cnv = py"mtlf"(times)
    return cnv
end