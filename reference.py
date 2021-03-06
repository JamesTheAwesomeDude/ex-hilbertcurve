from functools import wraps

def functionalize_mutator(f):
	@wraps(f)
	def g(x, *args, **kwargs):
		y = x.copy()
		f(y, *args, **kwargs)
		return y
	return g


def point(d, m, n=2):
	return undo_excess_work(gray_decode(transpose(d, m, n)), m, inverse=False)
	# return hilbertcurve.HilbertCurve(m, n).point_from_distance(d)


def dist(x, m):
	return untranspose(gray_encode(undo_excess_work(x, m, inverse=True), m), m)
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


@functionalize_mutator
def gray_decode(x):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L128-L132
	"""
	n = len(x)

	t = x[n-1] >> 1
	for i in range(n-1, 0, -1):
		x[i] ^= x[i-1]
	x[0] ^= t


def gray_encode(x, m=None):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L228-L238
	"""
	m = max(map(int.bit_length, x)) if m is None else m

	x = _gray_encode_a(x)
	t = _gray_encode_b(x, m=m)
	x = _gray_encode_c(x, t)

	return x


@functionalize_mutator
def _gray_encode_a(x):
	n = len(x)
	for i in range(1, n):
		x[i] ^= x[i-1]


def _gray_encode_b(x, m):
	n = len(x)
	t = 0
	q = 1 << (m - 1)
	while q > 1:
		p = q - 1
		if x[n-1] & q:
			t ^= p
		q >>= 1
	return t


@functionalize_mutator
def _gray_encode_c(x, t):
	n = len(x)
	for i in range(n):
		x[i] ^= t


def uew_atom(q, xi, x0):
	p = q - 1
	if not xi & q:
		# exchange low bits
		t = (x0 ^ xi) & p
		x0 ^= t
		xi ^= t
	else:
		# invert low bits
		x0 ^= p
	return xi, x0


def undo_excess_work(x, m=None, /, inverse=False):
	"""
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L134-L147
	https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L215-L226
	"""
	n = len(x)
	m = max(map(int.bit_length, x)) if m is None else m

	if not inverse:
		k_range = range(0, m, 1) # (0..(m-1))
	else:
		k_range = range(m-1, 0, -1) # ((m-1)..1)

	for q in map(lambda k: 1 << k, k_range):
		x = undo_excess_work_inner(x, q, inverse=inverse)
	return x


@functionalize_mutator
def undo_excess_work_inner(x, q, inverse):
	#print("undo_excess_work_inner(%s, %i, %s)" % (repr(x), q, repr(inverse)))
	n = len(x)
	if not inverse:
		i_range = range(n-1, -1, -1) # ((n-1)..0)
	else:
		i_range = range(0, n, 1) # (0..(n-1))

	for i in i_range:
		x[i], x[0] = uew_atom(q, x[i], x[0])

	#print(" => %s" % (repr(x),))
