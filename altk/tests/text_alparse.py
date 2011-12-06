import alparse as al

def test_replace_linear_mat():
    text = "aMat[0][1] = 12\nbMat[23][45] = x + z\n"
    matrixNames = ('aMat', 'bMat', 'cMat', 'dMat')
    replaced = al.replace_linear_mat(matrixNames, text)
    assert replaced == "A0, 1] = 12\nB[23, 45] = x + z\n"
