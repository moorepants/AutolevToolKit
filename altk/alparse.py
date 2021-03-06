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
    '''Sets the file location to the first line matching string. With reference
    to the beginning of the file.'''
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

def write_python(inFileStrings, cFileStrings, className, matrixNames, directory=None):
    '''Writes a basic Python class definition.

    '''

    if not directory == None:
        classFile = os.path.join(directory, className)
    else:
        classFile = className

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

    # open up the template file
    template = open(os.path.join(os.path.dirname(__file__), 'templates',
        'DynamicSystemTemplate.txt'), 'r')
    # grab the text from the template file
    data = template.read()
    template.close()

    # substitute for all the a tags
    stateNameLines, initCondLines = state_and_initial_lines(states)
    data = re.sub('<stateNames>', stateNameLines, data)
    data = re.sub('<initialConditions>', initCondLines, data)

    intOptsDict = variable_declarations_to_dictionary(intopts)
    intOptString = write_dictionary('intOpts', intOptsDict, indentation=4)
    data = re.sub('<intOpts>', intOptString, data)
    parDict = variable_declarations_to_dictionary(parameters)
    parameterString = write_dictionary('parameters', parDict, indentation=4)
    data = re.sub('<parameters>', parameterString, data)

    data = re.sub('<name>', className, data)

    data = re.sub('<inputNames>', input_lines(inputs)[0], data)
    data = re.sub('<inputs>', input_lines(inputs)[1], data)
    data = re.sub('<zeroInputs>', zero_inputs(inputs), data)

    stateNames, initialConditions = variables_values(states)
    oNames, oLines = output_lines(outputNames, outputs)
    data = re.sub('<outputNames>', oNames, data)

    data = re.sub('<outputs>', oLines, data)
    data = re.sub('<numZees>', zee_line(variables), data)
    inputNames, blah = variables_values(inputs)
    data = re.sub('<eom>', eom_lines(parDict, stateNames, inputNames, odefunc), data)
    data = re.sub('<constants>', constants_lines(constants), data)
    data = re.sub('<dependent>', dependentVarLines, data)

    data = re.sub('<kinematical>', self_dot_z(extract_kinematical(odefunc,
        stateNames)), data)

    constantNames, blah = variables_values(constants)

    data = re.sub('<extractParameters>',
            create_extract_parameter_lines(parDict.keys()), data)
    data = re.sub('<extractConstants>',
            create_extract_parameter_lines(constantNames), data)
    data = re.sub('<extractStates>', create_extract_state_lines(stateNames), data)
    data = re.sub('<linear>', self_dot_z(replace_linear_mat(matrixNames,
        indent(linear, 8))), data)

    # write the modified data to file
    outputfile = open(classFile + '.py', 'w')
    outputfile.write(data)
    outputfile.close()

def extract_kinematical(odefun, stateNames):
    kinematical = ''
    for line in odefun.splitlines():
        for state in stateNames:
            if line.startswith(state):
                kinematical += line + '\n'
    return indent(kinematical, 8)

def zero_inputs(inputs):
    """Returns a line which sets each input variable equal to zero."""
    zeroInputs = ''
    for inputVar in inputs.splitlines():
        zeroInputs += inputVar.split('=')[0].strip() + ' = 0.0\n'
    return indent(zeroInputs, 8)

def indent(text, indentation):
    """Adds the desired indentation to string of lines."""
    itext = re.sub(r'\n', r'\n' + ' ' * indentation, text)
    return ' ' * indentation + itext[:-indentation]

def replace_linear_mat(matrixNames, text):
    """Returns a string such that the linear matrix names are formatted for the
    python output."""

    for mat in zip(matrixNames, ('self.A', 'self.B', 'self.C', 'self.D')):
        text = re.sub(mat[0] + r'\[(\d*)\]\[(\d*)\]', mat[1] + r'[\1, \2]', text)
    return text

def first_line(string, numIndents):
    firstLine = ' ' * 4 * numIndents + string
    indent = len(firstLine)
    return firstLine, indent

def self_dot_z(string):
    '''Returns a string with z[x] changed to self.z[x].'''
    return re.sub('(z\[\d*\])', r'self.\1', string)

def write_list(varName, valList, indentation=0, oneLine=False):
    '''Returns a text string for a list declaration.

    Parameters
    ----------
    varName : string
        The name of the variable to store the list.
    valList : list
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
    for i, val in enumerate(valList):
        # if it is a string put quotes around it
        if type(val) == type('z'):
            val = "'" + val + "'"
        else:
            val = str(val)
        if i == 0:
            listString += val + afterVar
        elif i == len(valList) - 1:
            listString += ' '*indent + val + ']'
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
        and the value is a float.
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
            dictString += ' ' * indent + line + "}"
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

def create_extract_parameter_lines(parameterNames, indentSpaces=8):
    """Returns a string of lines which extract the parameters from the
    parameter dictionary.

    Parameters
    ----------
    parameters : list
        A list of the model parameters.

    Returns
    -------
    parameterLines : string
        A string of lines that will extract the parameters from the parameter
        dictionary.

    """
    indent = ' ' * indentSpaces

    # create the parameter declaration lines
    parameterLines = indent + '# declare the parameters\n'
    for k in parameterNames:
        parameterLines += indent + k + " = self.parameters['" + k + "']\n"

    return parameterLines

def create_extract_state_lines(stateNames, indentSpaces=8):
    """Returns a string of lines which extract the parameters from the
    parameter dictionary.

    Parameters
    ----------
    stateNames : list
        A list of the state names.

    Returns
    -------
    extractStateLines : string
        A string of lines that will extract the states from the state array.

    """
    indent = ' ' * indentSpaces

    # create the state declaration lines
    extractStateLines = indent + '# declare the states\n'
    for i, name in enumerate(stateNames):
        extractStateLines += indent + name + ' = x[' + str(i) + ']\n'

    return extractStateLines

def eom_lines(parameters, stateNames, inputNames, odefunc):
    """Returns the lines for the equations of motion section of the python
    file.

    Parameters
    ----------
    parameters : dictionary
        A dictionary of the model parameters.
    stateNames : list
        A list of the state names.
    inputNames : list
        A list of the input names.
    odefun : string
        A string containing the essential equations of motion of the system.

    Returns
    -------
    f : string
        A string containing the equations of motion setup for the function `f`
        in the DynamicSystem class.

    """

    indent = ' ' * 8

    # create the input declaration lines
    inputLines = indent + '# calculate and declare the inputs\n'
    inputLines += indent + 'u = self.inputs(t)\n'
    for i, name in enumerate(inputNames):
        inputLines += indent + name + ' = u[' + str(i) + ']\n'

    # create the equation of motion lines
    eom = odefunc.splitlines()
    eomLines = indent + '# calculate the derivatives of the states\n'
    # if there are zee's in the lines substute them with self.z[...]
    for line in eom:
        eomLines += indent + self_dot_z(line) + '\n'

    # create the derivatives lines
    derivativeLines = indent + '# store the results in f and return\n'
    derivativeLines += indent + 'f = zeros_like(x)\n'
    for i, name in enumerate(stateNames):
        derivativeLines += indent + 'f[' + str(i) + '] = ' + name +'p\n'

    f = inputLines + '\n'
    f += eomLines + '\n' + derivativeLines

    return f

def zee_line(variables):
    print "processing the zee number"
    # find the z variable declaration
    firstCharacters = [x[:2] for x in variables]
    try:
        index = firstCharacters.index('z[')
        numZees = re.sub('z\[(\d*)\]', r'\1', variables[index])
        return ' '*4 + 'z = zeros(' + numZees + ')'
    except:
        return '    # no zees here'

def output_lines(outputNames, outputs):
    print "processing the outputs"

    indent = ' ' * 8

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
        outputLines += self_dot_z(indent + line + '\n')

    # create the output declarations
    oDecLines = indent + '# store the results in y and return\n'
    oDecLines += indent + 'y = zeros(len(self.outputNames))\n'
    for i, name in enumerate(outputNames):
        oDecLines += indent + 'y[' + str(i) + '] = ' + name + '\n'

    outputLines = outputLines + '\n' + oDecLines

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

def variables_values(equationString):
    """Returns an ordered list of both the variable names and the values.

    """
    equationLines = equationString.splitlines()
    variableList = []
    valueList = []
    for line in equationLines:
        var, val = line.split('=')
        if var.startswith('z'):
            pass
        else:
            variableList.append(var.strip())
            try:
                valueList.append(float(val.strip()))
            except ValueError:
                valueList.append(val.strip())


    return variableList, valueList

def state_and_initial_lines(states):
    """Returns the

    """
    states = states.splitlines()
    stateNameList = []
    initCondList = []
    for state in states:
        var, val = state.split(' = ')
        stateNameList.append(var)
        initCondList.append(float(val))

    stateLines = write_list('stateNames', stateNameList, indentation=4)
    initLines = write_list('initialConditions', initCondList, indentation=4)
    return stateLines, initLines

def variable_declarations_to_dictionary(equationString):
    '''Write a multiline string of simple variable declarations as a
    dictionary.

    Parameters
    ----------
    equationString : string
        This string should be in the form "a = 2.0\nb = 3.0\n"

    Returns
    -------
    equationDict : dictionary

    '''
    equationList = equationString.splitlines()
    equationDict = {}
    for line in equationList:
        key, val = line.split('=')
        equationDict[key.strip()] = float(val.strip())

    return equationDict

def writeCxx(inFileStrings, cFileStrings, className):
    raise Exception

def alparsein(fileNameBase, code):
    """Parse the .in file from Autolev to grab all the lines that begin with
    the word 'Constant' or 'Initial Value'
    """

    print "cwd: ", os.getcwd()
    fp = open(fileNameBase + "Dynamics.in", "r")
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

# this function is also inside of state_and_initial_lines, they need to be
# merged
def state_name_list(states):
    stateNames = []
    for line in states.splitlines():
        state, val = line.split(' = ')
        stateNames.append(state)
    return stateNames

def equation_lines_to_dictionary(lines):
    '''Returns a dictionary such that the left hand side of the equations
    in lines is the keyword and the right hand side is the pair.

    lines should look like this: 'a = b\nc=  a + b + d\n' with an '=' as the
    separator.

    '''
    dictionary = {}
    for eq in lines.splitlines():
        var, expr = eq.split('=')
        dictionary[var.strip()] = expr.strip()
    return dictionary

def equation_dictionary_to_text(dictionary, code, sort=False):
    '''Returns a string with an equation on each line.'''
    keys = dictionary.keys()
    if sort:
        keys.sort()
    lines = ''
    if code == 'Python':
        ending = ''
    elif code == 'C':
        ending = ';'
    else:
        ending = ''
    for key in keys:
        lines += key + ' = ' + dictionary[key] + ending + '\n'
    return lines

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

    fp = open(fileNameBase + "Dynamics.c", "r")

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
            l += "\n"
            if foundSpecified:
                inputs += l
            else:
                odefunc += l

    # grab all the non zee equations out of the odefunc
    nonZees = []
    for eq in odefunc.splitlines():
        # store all of the equations that aren't zees
        if eq[0] != 'z':
            nonZees.append(eq)

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
                if l[:len(matrix) + 1] == matrix + '[':
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
            continue
        break

    linear = linearBeg + linear

    # write out the dependent variable equations
    dependentVarLines = ''
    for line in nonZees:
        var, rest = line.split(' = ')
        if var in outputNames:
            # add some indentation and replace the zees
            dependentVarLines += ' '*8 + self_dot_z(line) + '\n'

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

    # remove those stupid alTmp files in the directory
    try:
        os.system('rm ' + os.path.join(directory, 'alTmp.*'))
    except:
        pass

    inFileStrings = alparsein(fileNameBase, code)
    cFileStrings = alparsec(fileNameBase, code, linear,
                            state_name_list(inFileStrings[2]))

    if code == "Text":
        writeText(fileNameBase, className, inFileStrings, cFileStrings,
                      directory=directory)
    elif code == "C":
        writeC(inFileStrings, cFileStrings, className)
    elif code == "Python":
        write_python(inFileStrings, cFileStrings, className, linear, directory=directory)
    elif code == "C++":
        writeC++(inFileStrings, cFileStrings, className)

