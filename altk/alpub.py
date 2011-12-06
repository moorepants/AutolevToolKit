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
