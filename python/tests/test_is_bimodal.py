import is_bimodal

assert is_bimodal.is_histogram_bimodal([4, 1, 2, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 1]) == True 
assert is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 4, 1]) == True
assert is_bimodal.is_histogram_bimodal([1, 4, 1, 1, 4, 4, 1]) == False # should be True
assert is_bimodal.is_histogram_bimodal([1, 4, 2, 1, 1, 4, 1]) == True # should be True

print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4, 1]))
print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4, 4]))
print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4,]))
print(is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4]))
print(is_bimodal.van_der_eijk([1, 4, 1, 1, 1, 4]))
print(is_bimodal.van_der_eijk([4, 4, 1, 1, 4]))
print(is_bimodal.van_der_eijk([4, 4, 1, 1, 4, 4]))
print(is_bimodal.van_der_eijk([4, 4, 1, 1, 4, 4]))

assert is_bimodal.is_histogram_bimodal([4, 1, 2, 4, 1]) == True
