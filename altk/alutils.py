#! usr/bin/env python

def write_linearization(matrices, states, inputs, outputs, holonomic=None, filename=None):
    """Returns the text for jacobian calculations for the system matrices for
    Autolev.

    Parameters
    ----------
    matrices : tuple
        The desired variable names for the state, input, output and feed
        forward matrices.
    states : tuple
        The state variable names of the system.
    inputs : tuple
        The input variable names of the system.
    outputs : tuple
        The output variable names of the system.
    holonomic : tuple
        The name of the holonomic constraint equation and the independent
        coordinate.
    filename : string
        The path to a file where the resulting text will be saved.

    Returns
    -------
    text : string
        The input text for Autolev.

    """

    def derivative(matrix, i, j, row, col, holonomic, prime=False):
        if prime is True:
            prime = "'"
        else:
            prime = ''
        partial = (matrix + '[' + str(i + 1) + ', ' + str(j + 1) + '] = d(' +
                row + prime + ', ' + col + ')')
        if holonomic is not None:
            dependent = holonomic[0]
            constraint = holonomic[1]
            chainrule = (' + d(' + row + prime + ', ' + dependent + ') * d(' +
                constraint + ', ' + col + ') / d(' + constraint + ', ' + dependent
                + ')')
        else:
            chainrule = ''
        return partial + chainrule + '\n'

    text = ''
    # state matrix
    for i, row in enumerate(states):
        for j, col in enumerate(states):
            text += derivative(matrices[0], i, j, row, col, holonomic, prime=True)
        text += '\n'

    # input matrix
    for i, row in enumerate(states):
        for j, col in enumerate(inputs):
            text += derivative(matrices[1], i, j, row, col, holonomic, prime=True)
        text += '\n'

    # output matrix
    for i, row in enumerate(outputs):
        for j, col in enumerate(states):
            text += derivative(matrices[2], i, j, row, col, holonomic)
        text += '\n'

    # feed forward matrix
    for i, row in enumerate(outputs):
        for j, col in enumerate(inputs):
            text += derivative(matrices[3], i, j, row, col, holonomic)
        text += '\n'

    text += 'encode '
    for x in matrices[:-1]:
        text += x + ', '
    text += matrices[-1]

    if filename is not None:
        with open(filename, 'w') as f:
            f.write(text)

    return text
