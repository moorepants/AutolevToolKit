"""Parse C output file from Autolev code dynamics().

    Prints the following to stdout:

    Parameters with default values
    States variables listed with default initial conditions

    Section to evaluate constants
    Section to evaluate right hand side of state derivatives
    Section to evaluate output quantities

"""
import os
import re

def seekto(fp, string):
    '''Sets the file location to the first line matching string.'''
    # go back to the beginning of the file
    fp.seek(0)
    # go through line by line and stop when string is found
    for l in fp:
        if l.strip() == string:
            break

def writeText(fileNameBase, className, inFileStrings, cFileStrings,
        directory=None):
    '''Writes a plain text file with the model details.'''
    if not directory == None:
        classFile = os.path.join(directory, className)
    else:
        classFile = className
    intopts, parameters, states = inFileStrings
    variables, constants, odefunc, outputs, inputs, linear, outputNames, dependentVarLines = cFileStrings

    fp = open(classFile + ".txt", "w")

    fp.write("[Name]\n" + className + "\n\n")

    fp.write("[Integration Options]\n")
    fp.write(intopts + "\n")

    fp.write("[Parameters]\n")
    fp.write(parameters + "\n")

    fp.write("[States]\n")
    fp.write(states + "\n")

    fp.write("[Constants]\n")
    fp.write(constants + "\n")

    fp.write("[Inputs]\n")
    fp.write(inputs + "\n")

    fp.write("[Equations of Motion]\n")
    fp.write(odefunc + "\n")

    fp.write("[Dependent Variables]\n")
    fp.write(dependentVarLines + "\n")

    fp.write("[Output Names]\n")
    otnm = ''
    for name in outputNames:
        otnm += name + '\n'
    fp.write(otnm + "\n")

    fp.write("[Outputs]\n")
    fp.write(outputs + "\n")

    fp.write("[Linear]\n")
    fp.write(linear)

    print(fileNameBase + ".in and " + fileNameBase + ".c sucessfully" +
            " parsed.  Output code is in:\n" + fp.name)
    fp.close()

def writeC(inFileStrings, cFileStrings, className):
    raise Exception
    intopts, parameters, states = inFileStrings
    variables, constants, odefunc, outputs = cFileStrings

    fileNameBase += "_al"
    fp_header = open(fileNameBase + ".h", "w")
    fp_implementation = open(fileNameBase + ".c", "w")
    fp_driver = open(fileNameBase + "_main.c", "w")

    # Write the variables on one long line
    varstring = ""
    for v in variables:
        varstring += v + ", "
    varstring = varstring[:-2] + ";\n"

    # Write the header file
    fp_header.write(
        "#ifndef " + fileNameBase.upper() + "_H\n" +
        "#define " + fileNameBase.upper() + "_H\n\n" +
        "// All variables defined as globals with file scope\n" +
        "double " + varstring + "\n" +
        "// Function prototypes\n" +
        "int initConstants(void);\n" +
        "int eoms(double t, const double x[], double f[], void * p);\n" +
        "void outputs(void);\n" +
        "#endif")
    fp_header.close()

    # Write the implementation file
    indented_constants = ""
    indented_odefun = ""
    indented_outputs = ""
    for l in constants.splitlines(True):
        indented_constants += "  " + l
    for l in odefun.splitlines(True):
        indented_odefun += "  " + l
    for l in outputs.splitlines(True):
        indented_outputs += "  " + l

    fp_implementation.write(
        "#include <math.h>\n"
        "#include <gsl/gsl_odeiv.h>\n" +
        "#include \"" + fp_header.name() + "\"\n" +
        "int initConstants(void)\n{\n" + constants + "\n}" +
        "// initConstants()\n\n" +
        "int eoms(double t, const double x[], double f[], void * p)\n{\n")
        #for i in range(len(
        #
        #        "void outputs(void);\n" +
        #        "#endif")
    fp_implementation.close()

def writePython(inFileStrings, cFileStrings, className, directory=None):
    '''
    Writes a basic Python class definition.

    '''
    template = open('DynamicSystemTemplate.py', 'r')

    if not directory == None:
        classFile = os.path.join(directory, className)
    else:
        classFile = className

    outputfile = open(classFile + '.py', 'w')

    intopts, parameters, states = inFileStrings
    variables, constants, odefunc, outputs, inputs, linear, outputNames, dependentVarLines = cFileStrings
    #print "intopts:\n", intopts
    #print "parameters:\n", parameters
    #print "states:\n", states
    #print "variables:\n", variables
    #print "constants:\n", constants
    #print "odefunc:\n", odefunc
    #print "outputs:\n", outputs
    #print "inputs:\n", inputs
    #print "linear:\n", linear
    #print "outputNames:\n", outputNames

    state_lines(states)

    data = template.read()

    data = re.sub('<name>', className, data)
    data = re.sub('<stateNames>', state_lines(states)[0], data)
    data = re.sub('<initialConditions>', state_lines(states)[1], data)
    data = re.sub('<inputNames>', input_lines(inputs)[0], data)
    data = re.sub('<inputs>', input_lines(inputs)[1], data)
    data = re.sub('<outputNames>', output_lines(outputNames, outputs)[0], data)
    data = re.sub('<outputs>', output_lines(outputNames, outputs)[1], data)
    data = re.sub('<numZees>', zee_line(variables), data)
    data = re.sub('<eom>', eom_lines(odefunc), data)
    data = re.sub('<constants>', constants_lines(constants), data)
    data = re.sub('<dependent>', dependentVarLines, data)
    data = re.sub('<intOpts>', var_declarations_to_dictionary(intopts), data)
    data = re.sub('<parameters>', var_declarations_to_dictionary(parameters), data)
    outputfile.write(data)

    template.close()
    outputfile.close()

def first_line(string, numIndents):
    firstLine = ' '*4*numIndents + string
    indent = len(firstLine)
    return firstLine, indent

def self_dot_z(string):
    return re.sub('(z\[\d*\])', r'self.\1', string)

def write_list(varName, valList, indentation=0, oneLine=False):
    '''Returns a text string for a list declaration.

    Parameters
    ----------
    varName : string
        The name of the variable to store the list.
    valList : list of strings
        A list of the values to be stored in varName.
    indention : integer
        Number of space you want the string to be indented.
    oneLine : boolean
        Make this true if you want the list to be on one line.

    Returns
    -------
    listString : string
        A string formatted nicely to declare a list.

    '''
    listString = ' '*indentation + varName + ' = ['
    if oneLine:
        afterVar = ', '
        indent = 0
    else:
        indent = len(listString)
        afterVar = ',\n'
    for val in valList:
        if val == valList[0]:
            listString += val + afterVar
        elif val == valList[-1]:
            listString += ' '*indent + val + ']\n'
        else:
            listString += ' '*indent + val + afterVar
    return listString

def write_dictionary(varName, dictionary, indentation=0, oneLine=False):
    '''Returns a text string for a dictionary declaration.

    Parameters
    ----------
    varName : string
        The name of the variable to store the list.
    dictionary : dictionary
        A dictionary of the values to be stored in varName. The key is a string
        and the value is a number.
    indention : integer
        Number of space you want the string to be indented.
    oneLine : boolean
        Make this true if you want the list to be on one line.

    Returns
    -------
    dictString : string
        A string formatted nicely to declare a dictionary.

    '''
    dictString = ' '*indentation + varName + ' = {'
    if oneLine:
        afterVar = ', '
        indent = 0
    else:
        indent = len(dictString)
        afterVar = ',\n'
    keyList = dictionary.keys()
    keyList.sort()
    for key in keyList:
        line = "'" + key + "' : " + str(dictionary[key])
        if key == keyList[0]:
            dictString += line + afterVar
        elif key == keyList[-1]:
            dictString += ' ' * indent + line + "}\n"
        else:
            dictString += ' ' * indent + line + afterVar
    return dictString

def constants_lines(constants):
    print "processing constants"
    constants = constants.splitlines()
    constantsLines = ''
    constantList = []
    for line in constants:
        constantsLines += ' '*8 + self_dot_z(line) + '\n'
        if line[0] != 'z':
            var, trash = line.split(' = ')
            constantList.append(var)
    for cst in constantList:
        constantsLines = re.sub('(\W)(' + cst + ')(\W)',
                                r"\1self.parameters['\2']\3",
                                constantsLines)
    return constantsLines

def eom_lines(odefunc):
    print "processing the equations of motion"
    eom = odefunc.splitlines()
    eomLines = ''
    for line in eom:
        eomLines += ' '*8 + re.sub('(z\[\d*\])', r'self.\1', line) + '\n'
    return eomLines


def zee_line(variables):
    print "processing the zee number"
    # find the z variable declaration
    firstCharacters = [x[:2] for x in variables]
    print firstCharacters
    try:
        index = firstCharacters.index('z[')
        numZees = re.sub('z\[(\d*)\]', r'\1', variables[index])
        return ' '*4 + 'z = zeros(' + numZees + ')'
    except:
        return '    # no zees here'

def output_lines(outputNames, outputs):
    print "processing the outputs"
    outputNameLines, outputNameIndent = first_line('outputNames = [', 1)
    for name in outputNames:
        if name == outputNames[0]:
            outputNameLines += "'" + name + "',\n"
        elif name == outputNames[-1]:
            outputNameLines += outputNameIndent*' ' + "'" + name + "']"
        else:
            outputNameLines += outputNameIndent*' ' + "'" + name + "',\n"
    outputLines = ''
    for line in outputs.splitlines():
        outputLines += self_dot_z(' '*8 + line + '\n')
    return outputNameLines, outputLines

def input_lines(inputs):
    print "processing the inputs"
    inputs = inputs.splitlines()
    inputNameLines, inputNameIndent = first_line('inputNames = [', 1)
    for i, inpt in enumerate(inputs):
        var, expr = inpt.split(' = ')
        if inpt == inputs[0]:
            inputNameLines += "'" + var + "',\n"
            inputLines = ' '*8 + 'u[' + str(i) + '] = ' + expr + '\n'
        elif inpt == inputs[-1]:
            inputNameLines += inputNameIndent*' ' + "'" + var + "']"
            inputLines += ' '*8 + 'u[' + str(i) + '] = ' + expr
        else:
            inputNameLines += inputNameIndent*' ' + "'" + var + "',\n"
            inputLines += ' '*8 + 'u[' + str(i) + '] = ' + expr + '\n'
    return inputNameLines, inputLines

def state_lines(states):
    print "processing the states"
    states = states.splitlines()
    fsp = '    '
    stateLines = fsp + 'stateNames = ['
    stateIndent = len(stateLines)
    initLines = fsp + 'initialConditions = ['
    initIndent = len(initLines)
    for state in states:
        var, val = state.split(' = ')
        if state == states[0]:
            stateLines +=  "'" + var + "',\n"
            initLines += val + ',\n'
        elif state == states[-1]:
            stateLines +=  ' '*stateIndent + "'" + var + "']"
            initLines += ' '*initIndent + val + ']'
        else:
            stateLines +=  ' '*stateIndent + "'" + var + "',\n"
            initLines += ' '*initIndent + val + ',\n'
    return stateLines, initLines

def var_declarations_to_dictionary(mulLineEqs, varName, indentation=0, oneLine=False):
    '''Write a multiline string of simple variable declarations as a dictionary
    definiton for a Python file.

    Parameters
    ----------
    mulLineEqs : string
        This string should be in the form "a = 2.0\nb = 3.0\n"
    varName : string
        The name of the variable in the generated string.
    indentation : int
        Number of spaces of indentation for the entire string.
    oneLine : boolean
        Put the dictionary declaration all on one line.

    Returns
    -------
    varDecDict : dictionary

    '''
    eqs = mulLineEqs.splitlines()
    eqDict = {}
    for eq in eqs:
        key, val = eq.split(' = ')
        eqDict[key] = float(val)
    return write_dictionary(varName, eqDict,
                            indentation=indentation,
                            oneLine=oneLine)

def writeCxx(inFileStrings, cFileStrings, className):
    raise Exception

def alparsein(fileNameBase, code):
    """Parse the .in file from Autolev to grab all the lines that begin with
    the word 'Constant' or 'Initial Value'
    """

    print "cwd: ", os.getcwd()
    fp = open(fileNameBase + ".in", "r")
    for i in range(6):
        fp.next()

    intopts = ""
    parameters = ""
    states = ""

    for l in fp:
        l = l.strip().split()
        if l:
            if l[0] == "Constant":
                parameters += l[1] + " = " + l[4]
                if l[2] != "UNITS" and code == "Text":
                    parameters += ", " + l[2]
                if code == "C" or code == "C++":
                    parameters += ";"
                parameters += "\n"
            elif l[0] == "Initial" and l[1] == "Value":
                states += l[2] + " = " + l[5]
                if l[3] != "UNITS" and code == "Text":
                    states += ", " + l[3]
                if code == "C" or code == "C++":
                    states += ";"
                states += "\n"
            elif l[2] == 'TINITIAL':
                intopts += "ti = " + l[5]
                if code == "C" or code == "C++":
                    intopts += ";"
                elif code == "Text" and l[3] != "UNITS":
                    intopts += ", " + l[3]
                intopts += "\n"
            elif l[2] == 'TFINAL':
                intopts += "tf = " + l[5]
                if code == "C" or code == "C++":
                    intopts += ";"
                elif code == "Text" and l[3] != "UNITS":
                    intopts += ", " + l[3]
                intopts += "\n"
            elif l[2] == 'INTEGSTP':
                intopts += "ts = " + l[5]
                if code == "C" or code == "C++":
                    intopts += ";"
                elif code == "Text" and l[3] != "UNITS":
                    intopts += ", " + l[3]
                intopts += "\n"
            elif l[2] == 'ABSERR':
                intopts += "abserr = " + l[4]
                if code == "C" or code == "C++":
                    intopts += ";"
                intopts += "\n"
            elif l[2] == 'RELERR':
                intopts += "relerr = " + l[4]
                if code == "C" or code == "C++":
                    intopts += ";"
                intopts += "\n"
                break

    fp.close()
    return intopts, parameters, states

def state_name_list(states):
    stateNames = []
    for line in states.splitlines():
        state, val = line.split(' = ')
        stateNames.append(state)
    return stateNames

def alparsec(fileNameBase, code, linMat, stateNames):
    """Parse the .c file from Autolev to grab:
        1) list of variables that appear in all numerical calculations
        2) Evaluate constants section
        3) ode function
        4) output function
        5) specified inputs
        6) linear model ('A', 'B', 'C' and 'D' are reserver variable names)

        These 4 things are arranged in different ways, depending on the value
        of code and whether or not a class is to be automatically generated for
        C++ or Python code
    """

    fp = open(fileNameBase + ".c", "r")

    # For the Autolev C files I've examined, there are 20 lines of comments,
    # #include statements, and function forward declarations at the top.  The
    # following tosses these out the proverbial window.
    i = 0
    while i < 20:
        fp.next()
        i += 1

    # Loop to grab the statement that declares all the global variables,
    # assumes that they are declared as type 'double'
    variables = []
    for l in fp:
        l = l.strip().split()
        if l:
            if l[0] == "double":
                l = l[1].strip().split(',')

                if l[0] == "Pi":
                    l = l[1:]
                if l[0] == "DEGtoRAD":
                    l = l[1:]
                if l[0] == "RADtoDEG":
                    l = l[1:]

                # multi line statement
                while l[-1] == '':
                    l.pop(-1)
                    l += fp.next().strip().split(',')

                if l[-1][-1] == ';':
                    l[-1] = l[-1][:-1]

                # Get rid of the Encode[??]
                if l[-1][:6] == "Encode":
                    l.pop(-1)

            if l[0] == "/*" and l[2] == "MAIN" and l[4] == "*/":
                break
            variables += l

    # Seek to the line has the comment above the constants
    seekto(fp, "/* Evaluate constants */")

    # Get all the equations for the Evaluate constants
    constants = ""
    for l in fp:
        l = l.strip()
        if l:
            # Handle multi-line statements
            while l[-1] != ';':
                l += fp.next().strip()

            if code == "Text" or code == "Python":
                l = l[:-1]  # remove the semi-colon at end
            l += "\n"
            constants += l
        else:
            break

    # Seek to the line in the ode func that has the comment above the equations
    seekto(fp, "/* Update variables after integration step */" )
    while fp.next().strip() != '':
        continue

    # Get the equations in the right hand side of the odes
    odefunc = ""
    inputs = ""
    foundSpecified = False
    nonZees = []
    for l in fp:
        l = l.strip()
        if l == "/* Update derivative array prior to integration step */":
            break

        # Set to false when the next empty line is found after the specified
        # comment
        if foundSpecified and l == '':
            foundSpecified = False

        # The specified inputs are somewhere in the equations of motion
        sp1 = "/* Quantities to be specified */"
        sp2 = "/* Quantities which were specified */"

        if l == sp1 or l == sp2:
            foundSpecified = True
            pass # skip the line
        elif l == '':
            pass
        else:
            # Handle multi-line statements
            while l[-1] != ';':
                l += fp.next().strip()
            if code == "Text" or code == "Python":
                l = l[:-1]
            # store all of the equations that aren't zees
            if l[0] != 'z':
                nonZees.append(l)
            l += "\n"
            if foundSpecified:
                inputs += l
            else:
                odefunc += l
    print nonZees

    outputs = ""
    linear = ""
    outputNames = []
    ol = 'writef(Fptr['
    seekto(fp, "/* Write output to screen and to output file(s) */")
    for l in fp:
        l = l.strip()
        if l:
            # grab the output names
            if l[:len(ol)] == ol:
                outputNames += [x.strip() for x in l.split(',')[2:-1]]
            # skip the lines with writef
            if l[:6] == "writef":
                continue
            # Handle multi-line statements
            while l[-1] != ';':
                l += fp.next().strip()
            # removes the ';' at the end
            if code == "Text" or code == "Python":
                l = l[:-1]
            l += "\n"
            # if it is a matrix entry for the A, B, C, D matrices then put it
            # in the linear listing, else put it in the outputs section. This
            # section seems to typically only have the encoded matrices
            # anyways.
            linearLine = False
            for matrix in linMat:
                if l[:len(matrix)] == matrix:
                    linearLine = True
            if linearLine:
                linear += l
            else:
                outputs += l
            continue
        break

    # Seek to the first line of the output equations
    # The outputs seem to come before the zees associated with the encoded A,
    # B, C, D matrices
    seekto(fp, "/* Evaluate output quantities */")
    linearBeg = ''
    numOutputsFound = 0
    nonStateOutputs = []
    dependentVars = [x.split(' = ')[0] for x in nonZees]
    for name in outputNames:
        if name not in stateNames and name not in dependentVars:
            nonStateOutputs.append(name)
    print nonStateOutputs
    for l in fp:
        l = l.strip()
        if l:
            # Handle multi-line statements
            while l[-1] != ';':
                l += fp.next().strip()
            if code == "Text" or code == "Python":
                l = l[:-1]
            l += "\n"
            leftSide, rightSide = l.split(' = ')
            if numOutputsFound < len(nonStateOutputs):
                outputs += l
            else:
                linearBeg += l
            if leftSide in nonStateOutputs:
                numOutputsFound += 1
                print numOutputsFound
            continue
        break

    linear = linearBeg + linear

    # write out the dependent variable equations
    dependentVarLines = ''
    for line in nonZees:
        print line
        # add some indentation and replace the zees
        dependentVarLines += ' '*8 + self_dot_z(line) + '\n'

    print dependentVarLines

    stuff = (variables, constants, odefunc, outputs, inputs, linear,
            outputNames, dependentVarLines)

    return stuff

def alparse(fileNameBase, className, code="Text", directory=None,
            linear=('A','B','C','D')):
    """
        fileNameBase : string of the base input filename.  alparse() expects
        that fileNameBase.c and fileNameBase.in exist in the current working
        directory unless a directory is specified.

        className : Name of system.  Used to name classes in
        Psuedo-Code/C++/Python code, used to name struct in C code.  Output
        code is written to a file of title className. It should be in camel
        case.*

        code : valid choices are "Text", "Python", "C" or "C++"

        directory : Optional path to the directory in which fileNameBase.c and
        fileNameBase.in exist. If 'None', alparse assumes files are in current
        working directory.

        linear : These should reflect what the linear A, B, C, D matrices in
        the autolev file are named, change them if you use variables other than
        the default.

    """
    if not directory == None:
        fileNameBase = os.path.join(directory, fileNameBase)

    # remove those stupid alTmp files
    os.system('rm ' + os.path.join(directory, 'alTmp.*'))

    inFileStrings = alparsein(fileNameBase, code)
    cFileStrings = alparsec(fileNameBase, code, linear,
                            state_name_list(inFileStrings[2]))

    if code == "Text":
        writeText(fileNameBase, className, inFileStrings, cFileStrings,
                      directory=directory)
    elif code == "C":
        writeC(inFileStrings, cFileStrings, className)
    elif code == "Python":
        writePython(inFileStrings, cFileStrings, className, directory=directory)
    elif code == "C++":
        writeC++(inFileStrings, cFileStrings, className)

