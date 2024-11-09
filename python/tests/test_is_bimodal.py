import is_bimodal

def test_simple():
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


def test_hist(hist_file):
    hists = []
    with open(hist_file) as f:
        for line in f:
            fs = [x.strip() for x in line.split(',') if x.strip() ]
            expected = bool(fs[0])
            hist = [ int(x) for x in fs[1:] ]
            estimate = is_bimodal.is_histogram_bimodal(hist)
            hists.append((expected, estimate, hist))

    print(hists)

if __name__ == '__main__':
    test_simple()
    test_hist("../../rust/test/hist.txt")
