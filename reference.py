


def point(d, m, n=2):
	return undo_excess_work(gray_decode(transpose(d, m, n)), m)
	# return hilbertcurve.HilbertCurve(m, n).point_from_distance(d)


def dist(x, m):
	return untranspose(gray_encode(inverse_undo_excess_work(x, m), m), m)
	# return hilbertcurve.HilbertCurve(m, n).distance_from_point(x)


def transpose(d, m, n):
	d_str = format(d, 'b').zfill(m * n)
	return [int(d_str[i::n], 2) for i in range(n)]
	# return hilbertcurve.HilbertCurve(m, n)._hilbert_integer_to_transpose(d)


def untranspose(x, m):
	n = len(x)
	x_str = [format(x[i], 'b').zfill(m) for i in range(n)]
	return int(''.join([y[i] for i in range(m) for y in x_str]), 2)
	# return hilbertcurve.HilbertCurve(m, n)._transpose_to_hilbert_integer(x)


def gray_decode(x):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L128-L132
	"""
	x = x.copy()
	n = len(x)

	t = x[n-1] >> 1
	for i in range(n-1, 0, -1):
		x[i] ^= x[i-1]
	x[0] ^= t

	return x


def gray_encode(x, m=None):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L228-L238
	"""
	x = x.copy()
	n = len(x)
	if m is None:
		m = max(map(int.bit_length, x))

	for i in range(1, n):
		x[i] ^= x[i-1]
	t = 0
	q = 1 << (m - 1)
	while q > 1:
		if x[n-1] & q:
			t ^= q - 1
		q >>= 1
	for i in range(n):
		x[i] ^= t

	return x

def undo_excess_work(x, m=None, n=None):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L134-L147
	"""
	x = x.copy()
	if n is None:
		n = len(x)
	if m is None:
		m = max(map(int.bit_length, x))

	z = 2 << (m - 1)
	q = 2
	while q != z:
		p = q - 1
		for i in range(n-1, -1, -1):
			if x[i] & q:
				# invert
				x[0] ^= p
			else:
				# exchange
				t = (x[0] ^ x[i]) & p
				x[0] ^= t
				x[i] ^= t
		q <<= 1

	return x


def inverse_undo_excess_work(x, m=None, n=None):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L215-L226
	"""
	x = x.copy()
	if n is None:
		n = len(x)
	if m is None:
		m = max(map(int.bit_length, x))

	q = 1 << (m - 1)
	while q > 1:
		p = q - 1
		for i in range(n):
			if not x[i] & q:
				# exchange
				t = (x[0] ^ x[i]) & p
				x[0] ^= t
				x[i] ^= t
			else:
				# invert
				x[0] ^= p
		q >>= 1

	return x
