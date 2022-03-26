#!/usr/bin/env elixir

defmodule HilbertCurve do
	use Bitwise

	# Any function without a `@doc` is private.
	# Any private function that's declared with `def` is broken
	# or in need of testing / fixing-up.
	# Once every function is (either `defp` or `@doc`'d),
	# the module will be ready to use.

	@doc """
	Returns the co-ordinate list of an `ndims`-dimensional point `d` units
	along an `order`-degree hilbert curve approximation
	"""
	def point(d, order, ndims \\ 2) do
		# side_length = 2 ** order
		# volume = side_length ** ndims
		d
		|> transpose(order, ndims)
		|> gray_decode()
		|> uew(order, :normal)
	end

	@doc """
	Given a co-ordinate list `x` referring to a point, returns its distance
	along an `order`-degree hilbert curve approximation
	"""
	def dist(x, order) when is_list(x) do
		x
		|> uew(order, :invert)
		|> gray_encode(order)
		|> untranspose(order)
	end

	def gray_encode(x, order) do
		# https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L228-L238
		x = gray_encode_a(x)
		t = gray_encode_b(x, order)
		x = gray_encode_c(x, t)

		x
	end

	def gray_encode_a(x) do
		# for i in range(1, self.n):
		#   point[i] ^= point[i-1]
		Enum.map(Stream.zip(x, [0 | x]), fn {a, b} -> Bitwise.bxor(a, b) end)
	end

	def gray_encode_b(x, order) do
		# q = 1 << (order - 1)
		# while q > 1:
		#   if point[self.n-1] & q:
		#     t ^= q - 1
		#   q >>= 1
		x_last = List.last(x)
		max_k = order - 1
		Enum.reduce(
			(max_k..2),
			(0),
			fn cur, acc ->
				k = cur
				t = (acc)
				q = 1 <<< k
				p = q - 1
				if Bitwise.band(x_last, q) != 0 do
					(Bitwise.bxor(t, p))
				else
					(t)
				end
			end
		)
	end

	def gray_encode_c(x, t) do
		# for i in range(self.n):
		#   point[i] ^= t
		Enum.map(x, fn x_i -> Bitwise.bxor(x_i, t) end)
	end

	defp gray_decode(x) do
		# https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L128-L132
		t = List.last(x) >>> 1
		x0 = Bitwise.bxor(hd(x), t)
		[x0 | Enum.map(Stream.zip(x, tl x), fn {a, b} -> Bitwise.bxor(a, b) end)]
	end

	def uew_atom(q, {x_i, x_0}) do
		p = q - 1
		if Bitwise.band(x_i, q) == 0 do
			# Exchange low bits of x[i] and x[0]
			exchange_bits({x_i, x_0}, p)
		else
			# Just invert low bits of x[0]
			{x_i, invert_bits(x_0, p)}
		end
	end

	def uew(x, order, mode) do
		invert = case mode do :normal -> false; :invert -> true; end
		k_max = order - 1
		k_range = conditional_reverse(invert, (1..k_max))
		# for( q = 2 ; q != (2 << (order - 1)) ; q <<= 1 )
		# for( q = (1 << (order - 1)) ; q > 1 ; q >>= 1 )
		{x} = Enum.reduce(
			k_range,
			{x},
			fn cur, acc ->
				k = cur
				{x} = acc

				q = 1 <<< k

				x = scan_unzip(
					conditional_reverse(invert, x),
					{nil, hd(x)},
					fn cur, acc ->
						x_i = cur
						{_, x_0} = acc
						uew_atom(q, {x_i, x_0})
					end
				)
				|> (fn {l, _} -> conditional_reverse(invert, x) end).()

				{x}
			end
		)

		x # return
	end

	defp invert_bits(i, mask \\ -1) do
		Bitwise.bxor(i, mask);
	end

	defp exchange_bits({a, b}, mask) do
		c = Bitwise.bxor(a, b) &&& mask
		{Bitwise.bxor(a, c), Bitwise.bxor(b, c)}
	end

	defp conditional_reverse(condition, l) do
		if(condition, do: l, else: Enum.reverse(l))
	end

	defp transpose(i, order, ndims) do
		# https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L85-L97
		#TODO: figure out what this function does
		i
		|> to_bits(order * ndims)
		|> Stream.chunk_every(ndims)
		|> Stream.zip()
		|> Stream.map(&Tuple.to_list/1)
		|> Enum.map(&from_bits/1)
	end

	defp untranspose(x, order) do
		# https://github.com/galtay/hilbertcurve/blob/v2.0.5/hilbertcurve/hilbertcurve.py#L100-L112
		#TODO: figure out what this function does
		x
		|> Stream.map(fn i -> to_bits(i, order) end)
		|> Stream.zip()
		|> Enum.map(&Tuple.to_list/1)
		|> List.flatten()
		|> from_bits()
	end

	defp scan_unzip(enumerable, initial_acc, f) do
		Stream.scan(enumerable, initial_acc, f)
		# |> Enum.to_list() |> IO.inspect()
		|> Enum.unzip()
	end

	defp to_bits(i, width) do
		i
		|> to_bitstring(width)
		|> String.to_charlist()
		|> Enum.map(&(&1 - hd '0'))
	end

	defp to_bitstring(i, width \\ 0) do
		Integer.to_string(i, 2)
		|> String.pad_leading(width, "0")
	end

	# defp from_bits(bs) when is_bitstring(bs) do
	#	 padding_size = 7 - rem(:erlang.bit_size(bs) + 7, 8)
	#	 padding = <<0::size(padding_size)>>
	#	 :binary.decode_unsigned(<< padding, bs :: bitstring >>)
	# end

	defp from_bits(l) when is_list(l) do
		Enum.map(l, &(&1 + hd '0'))
		|> List.to_string()
		|> from_bitstring()
	end

	defp from_bitstring(s) do
		Integer.parse(s, 2)
		|> (fn {i, ""} -> i end).()
	end

	def tests do

		if [5, 10, 20] == point(7865, 5, 3) do
			IO.puts("HilbertCurve.point OK")
		else
			IO.puts("HilbertCurve.point BROKEN")

			# transpose function IS WORKING
			if [10, 14, 27] == transpose(7865, 5, 3) do
				IO.puts("HilbertCurve.transpose OK")
			else
				IO.puts("HilbertCurve.transpose BROKEN")
			end

			# gray_decode function IS WORKING
			if [7, 4, 21] == gray_decode([10, 14, 27]) do
				IO.puts("HilbertCurve.gray_decode OK")
			else
				IO.puts("HilbertCurve.gray_decode BROKEN")
			end

			# uew function IS BROKEN
			if [5, 10, 20] == uew([7, 4, 21], 5, :normal) do
				IO.puts("HilbertCurve.uew :normal OK")
			else
				IO.puts("HilbertCurve.uew :normal BROKEN")
			end

		end

		if 7865 == dist([5, 10, 20], 5) do
			IO.puts("HilbertCurve.dist OK")
		else
			IO.puts("HilbertCurve.dist BROKEN")
	
			# untranspose function IS WORKING
			if 7865 == untranspose([10, 14, 27], 5) do
				IO.puts("HilbertCurve.untranspose OK")
			else
				IO.puts("HilbertCurve.untranspose BROKEN")
			end

			# gray_encode function IS BROKEN
			if [10, 14, 27] == gray_encode([7, 4, 21], 5) do
				IO.puts("HilbertCurve.gray_encode OK")
			else
				IO.puts("HilbertCurve.gray_encode BROKEN")
			end

			# inverse uew function IS BROKEN
			if [7, 4, 21] == uew([5, 10, 20], 5, :invert) do
				IO.puts("HilbertCurve.uew :invert OK")
			else
				IO.puts("HilbertCurve.uew :invert BROKEN")
			end

		end

	end

end
