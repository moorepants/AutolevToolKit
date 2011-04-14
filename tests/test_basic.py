import filecmp
import alparse as alp
import os

def test_test1_al():
    alp.alparse("test1_al", "test1_al")
    assert filecmp.cmp('test1_al.txt', 'test1_al_desired.txt')
    #os.remove('test1_al.txt')

def test_write_list():
    oneLine = '  x = [1, 2, 3, 4]'
    mulLine = '  x = [1,\n       2,\n       3,\n       4]'
    oneResult = alp.write_list('x', [1, 2, 3, 4], indentation=2, oneLine=True)
    mulResult = alp.write_list('x', [1, 2, 3, 4], indentation=2)
    print oneLine, oneResult
    print mulLine, mulResult
    assert oneResult == oneLine
    assert mulResult == mulLine

def test_write_dictionary():
    oneLine = "  x = {'a' : 1.0, 'b' : 2.0, 'c' : 3.0}"
    mulLine = "  x = {'a' : 1.0,\n       'b' : 2.0,\n       'c' : 3.0}"
    oneResult = alp.write_dictionary('x', {'a':1.0,'b':2.0,'c':3.0},
                                     indentation=2, oneLine=True)
    mulResult = alp.write_dictionary('x', {'a':1.0,'b':2.0,'c':3.0},
                                     indentation=2)
    print oneLine, oneResult
    print mulLine, mulResult
    assert oneResult == oneLine
    assert mulResult == mulLine

def test_equation_lines_to_dictionary():
    lines = 'a = b\nc = a + b + d\nb=2*3*b/5\n'
    dictionary = {'a':'b','c':'a + b + d','b':'2*3*b/5'}
    result = alp.equation_lines_to_dictionary(lines)
    for key, val in result.items():
        try:
            assert dictionary[key] == val
        except KeyError:
            assert 1 == 2
