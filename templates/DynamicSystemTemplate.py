import os
from numpy import zeros, zeros_like
from numpy import sin, cos, tan
from alparse.DynamicSystem import DynamicSystem, LinearDynamicSystem

class <name>(DynamicSystem):
    """
<description>

    """

    # model name
    name = '<name>'

    filename = ''.join(name.split())
    directory = os.path.join('..', 'models', filename)

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
    t = intOpts['ti']

    def __init__(self):
        '''Just sets the constants.'''
        self.constants()

    def constants(self):
        '''Sets the zees that are constant.'''
<extractParameters>

<constants>
    def f(self, x, t):
        '''Returns the time derivative of the state vector.

        Parameters
        ----------
        x : ndarray, shape(n,)
            The state vector at this time.
        t : float
            Time.

        Returns
        -------
        f : ndarray, shape(n,)
            The time derivative of the state vector.

        '''
<extractParameters>

<extractConstants>

<extractStates>

<eom>
        return f

    def inputs(self, t):
        '''Returns the inputs to the system.

        Parameters
        ----------
        t : float
            Time.

        Returns
        -------
        u : ndarray, shape(m,)
            The input array as a function of time.

        '''
        T = t # this is hack because autolev likes to capitlize everything
        # initialize the u vector
        u = zeros(len(self.inputNames))
        # calculate the inputs
<inputs>
        return u

    def outputs(self, x):
        '''Returns the outputs of the system.

        Parameters:
        -----------
        x : ndarray, shape(n,)
            Current state

        Returns:
        --------
        y : ndarray, shape(m,)
            y(t)

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''
<extractParameters>

<extractConstants>

<extractStates>
        # these are dependent variables that may be needed for the main
        # calculations
<dependent>
        # calculate the outputs
<outputs>

        return y

class Linear<name>(<name>, LinearDynamicSystem):

    name = "Linear<name>"

    def inputs(self, t):
        '''Returns the inputs to the system.

        Parameters
        ----------
        t : float
            Time.

        Returns
        -------
        u : ndarray, shape(m,)
            The input array as a function of time.

        '''
        T = t # this is hack because autolev likes to capitlize everything
        # initialize the u vector
        u = zeros(len(self.inputNames))
        # calculate the inputs
<inputs>
        return u

    def f(self, x, t):
        '''Returns the derivative of the states'''

        u = self.inputs(t)

        xd = dot(self.A, x) + dot(self.B, u)

        return xd

    def outputs(self, x):

        u = self.inputs(t)

        y = dot(self.C, x) + dot(self.B, u)

        return y

    def linear(self, x):
        """Calculates the state, input, output and feedforward  matrices for the
        system linearized about the provided equilibrium point.

        Parameters
        ----------
        x : ndarray, shape(n,)
            The point at which to linearize the system about. The order of the
            values corresponds to the system states.

        """

        self.equilibriumPoint = x
        # sets the zees for the equilbrium points
        nonlin = <name>()
        nonlin.f(x, 0.)
        self.z = nonlin.z

<extractParameters>

<extractConstants>

<extractStates>

<dependent>

        # turns out that the qdots may be in the linear equations below, so you
        # really need the kinematical differential equations here so they will
        # be available

<kinematical>

        # the A needs to be self.A and needs to be initialize with zeros

        self.A = zeros((len(self.stateNames), len(self.stateNames)))
        self.B = zeros((len(self.stateNames), len(self.inputNames)))
        self.C = zeros((len(self.outputNames), len(self.stateNames)))
        self.D = zeros((len(self.outputNames), len(self.inputNames)))

        # also the inputs show up in the following equations, not sure what
        # that means

<zeroInputs>

<linear>
