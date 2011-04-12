import filecmp
import alparse as alp
import os

def test_test1_al():
    alp.alparse("test1_al", "test1_al")
    assert filecmp.cmp('test1_al.txt', 'test1_al_desired.txt')
    #os.remove('test1_al.txt')

def test_write_list():
    oneLine = '  x = [1, 2, 3, 4]\n'
    mulLine = '  x = [1,\n       2,\n       3,\n       4]\n'
    assert alp.write_list('x', ['1', '2', '3', '4'], indentation=2, oneLine=True) == oneLine
    assert alp.write_list('x', ['1', '2', '3', '4'], indentation=2) == mulLine

def test_write_dictionary():
    oneLine = "  x = {'a' : 1.0, 'b' : 2.0, 'c' : 3.0}\n"
    mulLine = "  x = {'a' : 1.0,\n       'b' : 2.0,\n       'c' : 3.0}\n"
    assert alp.write_dictionary('x', {'a':1.0,'b':2.0,'c':3.0},
                                indentation=2, oneLine=True) == oneLine
    assert alp.write_dictionary('x', {'a':1.0,'b':2.0,'c':3.0},
                                indentation=2) == mulLine
