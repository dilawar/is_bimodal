import math

import is_bimodal

assert is_bimodal.is_histogram_bimodal([4, 1, 2, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([1, 4, 1, 1, 4, 4, 1]) == False  # should be True
assert is_bimodal.is_histogram_bimodal([1, 4, 2, 1, 1, 4, 1]) == True  # should be True

# The worked example from van der Eijk (2001, p. 331).
assert math.isclose(
    is_bimodal.van_der_eijk([30, 40, 210, 130, 530, 50, 10]), 0.6113333, abs_tol=1e-6
)
assert is_bimodal.is_histogram_bimodal([30, 40, 210, 130, 530, 50, 10]) == False

# A perfectly uniform histogram has A == 0.
assert is_bimodal.van_der_eijk([5, 5, 5, 5]) == 0.0

# Peaks in the middle of the histogram are detected.
assert is_bimodal.is_histogram_bimodal([1, 4, 1, 1, 1, 4]) == True
assert is_bimodal.is_histogram_bimodal([4, 4, 1, 1, 4]) == True

# Bimodal histograms with peaks at the extremes.
assert is_bimodal.is_histogram_bimodal([10, 10, 0, 0, 0, 10, 10]) == True
assert is_bimodal.is_histogram_bimodal([10, 10, 0, 0, 0, 0, 10]) == True

# A single dominant peak is unimodal.
assert is_bimodal.is_histogram_bimodal([1, 1, 1, 1, 1, 1, 1000]) == False

print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4, 1]))
print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4, 4]))
print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4]))
print(is_bimodal.van_der_eijk([1, 4, 1, 1, 1, 4]))
print(is_bimodal.van_der_eijk([4, 4, 1, 1, 4]))
print(is_bimodal.van_der_eijk([4, 4, 1, 1, 4, 4]))

assert is_bimodal.is_histogram_bimodal([4, 1, 2, 4, 1]) == True

# Bimodality coefficient (Sarle/SAS).
assert math.isclose(
    is_bimodal.bimodality_coefficient([1, 1, 1, 1, 1]), 0.5882352941176471, abs_tol=1e-12
)
assert math.isclose(
    is_bimodal.bimodality_coefficient([10, 0, 0, 0, 10]), 1.0, abs_tol=1e-9
)
assert is_bimodal.bimodality_coefficient([0, 0, 10, 0, 0]) == 0.0

# Invalid input is rejected with a ValueError.
for bad in ([], [1, 2], [0, 0, 0]):
    for fn in (
        is_bimodal.van_der_eijk,
        is_bimodal.is_histogram_bimodal,
        is_bimodal.bimodality_coefficient,
    ):
        try:
            fn(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for {fn.__name__}({bad})")
