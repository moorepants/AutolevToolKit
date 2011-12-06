#! usr/bin/env python

"""This is crude parser for Autolev output files (typical extension is .all).
It converts the output lines to LaTeX compatible strings."""

import re

def all_parse(inputFile, outputFile, subMag=False):
    """Goes through the input file and writes new file with the output lines
    on a single line which is formatted for LaTeX."""

    allFile = open(inputFile, 'r')
    outputText = ''
    isOutputLine = False
    for line in allFile:
        # match lines that start with (\d and empty lines
        if re.match(r'\(\d', line.strip()) or line.strip() == '':
            # if the previous line was an output line then convert to LaTeX
            if isOutputLine is True:
                outputText += to_latex(previousLine, subMag=subMag) + '\r\n' + line
            # otherwise just copy the line
            else:
                outputText += line
            isOutputLine = False
        else:
            # output lines can be more than one line long and the first line
            # always starts with a '->'
            if line.startswith('->'):
                # start a line to add all the following lines to
                previousLine = line.strip()
                isOutputLine = True
            elif isOutputLine is True:
                previousLine += line.strip()
    allFile.close()

    # write the new file
    with open(outputFile, 'w') as f:
        f.write(outputText)

def to_latex(text, subMag=False):
    """This changes very specific Autolev output code to LaTeX."""
    # change COS(q1) and COS(q1) to c_1, s_1
    text = re.sub(r'SIN\(\w(\d)\)', r's_\1', text)
    text = re.sub(r'COS\(\w(\d)\)', r'c_\1', text)
    text = re.sub(r'SIN', r'sin', text)
    text = re.sub(r'COS', r'cos', text)
    # change unit vectors to hats, e2> to \hat{e}_2
    text = re.sub(r'(\w)(\d)>', r'\hat{\1}_\2', text)
    # dots instead of primes, u4' to \dot{u}_4
    text = re.sub(r"u(\d)'", r'\dot{u}_\1', text)
    # subscript variables
    text = re.sub(r'([a-zA-Z])(\d)', r'\1_\2', text)
    # wheel radii subscripts
    text = re.sub(r'r([RF])', r'r_\1', text)
    # w_a_n> to ^N\omega^A
    def low(match):
        return '^' + match.group(2).upper() + '\omega^' + match.group(1).upper()
    text = re.sub(r'W_([a-z])_([a-z])>', low, text)
    # remove the multiplication
    text = re.sub(r'\*', r'', text)
    if subMag is True:
        text = text.replace('(c_4^2c_5^2+(s_4s_7-s_5c_4c_7)^2)^0.5', 'm')
        #text = re.sub(r'\(c_4\^2c_5\^2+\(s_4s_7-s_5c_4c_7\)\^2\)\^\{0\.5\}', 'm', text)
    # 
    text = re.sub(r'\^(\d).(\d)', r'^{\1.\2}', text)

    return text

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
