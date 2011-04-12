import filecmp
import alparse as alp
import os

def test_test1_al():
    alp.alparse("test1_al", "test1_al")
    assert filecmp.cmp('test1_al.txt', 'test1_al_desired.txt')
    #os.remove('test1_al.txt')

def test_write_list():
    oneLine = '  x = [1, 2, 3, 4]\n'
    multLine = '  x = [1,\n       2,\n       3,\n       4]\n'
    assert alp.write_list('x', ['1', '2', '3', '4'], indentation=2, oneLine=True) == oneLine
    assert alp.write_list('x', ['1', '2', '3', '4'], indentation=2) == multLine
