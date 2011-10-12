import os
from numpy import zeros, zeros_like
from numpy import sin, cos, tan
from alparse.DynamicSystem import DynamicSystem

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
