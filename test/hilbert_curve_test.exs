defmodule HilbertCurveTest do
	use ExUnit.Case
	doctest HilbertCurve

	test "stock example" do
		# In 3-dimensional space,
		# the 7856th point along a
		# hilbert curve with 5 bits of precision
		# is x=5 y=10 z=20
		assert HilbertCurve.point(7865, 5, 3) == [5, 10, 20]
		assert HilbertCurve.dist([5, 10, 20], 5) == 7865
	end
end
