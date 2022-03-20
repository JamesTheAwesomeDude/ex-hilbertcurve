defmodule HilbertCurveTest do
  use ExUnit.Case
  doctest HilbertCurve

  test "stock example" do
    assert HilbertCurve.point(7865, 5, 3) == [5, 10, 20]
    assert HilbertCurve.point([5, 10, 20], 5) == 7865
  end
end
