from numpy import zeros
from DynamicSystem import DynamicSystem

class <name>(DynamicSystem):
    """
<description>

    """

    # model name
    name = '<name>'

    filename = ''.join(name.split())
    directory = 'models/' + filename + '/'

    # numerical integration parameters
<intOpts>

    # parameter names and their values
<parameters>

    # state names
<stateNames>

    # sets the initial conditions of the states
<initialConditions>

    # input names
<inputNames>

    # output names
<outputNames>

    # initialize state vector
    x = zeros(len(stateNames))

    # initialize output vector
    y = zeros(len(outputNames))

    # initialize input vector
    u = zeros(len(inputNames))

    # initializes the zees
<numZees>

    # intialize the time
    t = 0.0

    def __init__(self):
        '''Just sets the constants.'''
        self.constants()

    def constants(self):
        '''Sets the zees that are constant.'''
        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

<constants>
    def f(self, x, t):
        '''Returns the derivative of the states.

        Parameters:
        -----------
        x : ndarray
            State vector
        t : ndarray
            Time

        Returns:
        --------
        f : ndarray
            dx/dt

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''

        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        for i, name in enumerate(self.stateNames):
            exec(name + ' = ' + 'x[' + str(i) + ']')

        # calculates inputs
        u = self.inputs(t)
        for i, name in enumerate(self.inputNames):
            exec(name + ' = ' + 'self.u[' + str(i) + ']')

<eom>
        # plug in the derivatives for returning
        f = zeros(len(stateNames))
        for i, name in enumerate(self.stateNames):
            exec('f[' + str(i) + '] = ' + name + 'p')

        return f

    def inputs(self, t):
        '''Returns the inputs to the system.

        Parameters:
        -----------
        t : ndarray
            Time

        Returns:
        --------
        u : ndarray
            u(t)

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''
        t = T
        u = zeros(len(self.inputNames))
<inputs>
        return u

    def outputs(self, x):
        '''Returns the outputs of the system.

        Parameters:
        -----------
        x : ndarray
            Current state

        Returns:
        --------
        y : ndarray
            y(t)

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''
        # defines the parameters locally from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        for i, name in enumerate(self.stateNames):
            exec(name + ' = ' + 'x[' + str(i) + ']')

        # calculate the outputs
<outputs>
        # plug in the derivatives for returning
        y = zeros(len(self.outputNames))
        for i, name in enumerate(self.outputNames):
            exec('y[' + str(i) + '] = ' + name)
        return y
