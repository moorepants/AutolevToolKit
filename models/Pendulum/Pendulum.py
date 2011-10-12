import os
from numpy import zeros, zeros_like
from numpy import sin, cos, tan
from alparse.DynamicSystem import DynamicSystem

class Pendulum(DynamicSystem):
    """
<description>

    """

    # model name
    name = 'Pendulum'

    filename = ''.join(name.split())
    directory = os.path.join('..', 'models', filename)

    # numerical integration parameters
    intOpts = {'abserr' : 1e-08,
               'relerr' : 1e-07,
               'tf' : 1.0,
               'ti' : 0.0,
               'ts' : 0.1}

    # parameter names and their values
    parameters = {'g' : 9.81,
                  'i' : 0.5,
                  'l' : 2.0,
                  'm' : 4.0}

    # state names
    stateNames = ['omega',
                  'theta']

    # sets the initial conditions of the states
    initialConditions = [0.0,
                         0.0]

    # input names
    inputNames = ['torque',
                  'force']

    # output names
    outputNames = ['omega',
                   'theta',
                   'k',
                   'p',
                   'longoutput']

    # initialize state vector
    x = zeros(len(stateNames))

    # initialize output vector
    y = zeros(len(outputNames))

    # initialize input vector
    u = zeros(len(inputNames))

    # initializes the zees
    z = zeros(18)

    # intialize the time
    t = intOpts['ti']

    def __init__(self):
        '''Just sets the constants.'''
        self.constants()

    def constants(self):
        '''Sets the zees that are constant.'''
        # declare the parameters
        i = self.parameters['i']
        m = self.parameters['m']
        l = self.parameters['l']
        g = self.parameters['g']


        self.z[8] = g*m
        self.z[15] = g*l*m

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
        # declare the parameters
        i = self.parameters['i']
        m = self.parameters['m']
        l = self.parameters['l']
        g = self.parameters['g']


        # declare the parameters


        # declare the states
        omega = x[0]
        theta = x[1]


        # calculate and declare the inputs
        u = self.inputs(t)
        torque = u[0]
        force = u[1]

        # calculate the derivatives of the states
        thetap = omega
        self.z[1] = cos(theta)
        self.z[2] = sin(theta)
        self.z[3] = pow(self.z[1],2) + pow(self.z[2],2)
        self.z[4] = l*self.z[3]
        self.z[10] = i*self.z[3]
        self.z[11] = self.z[3]*self.z[10] + 0.25*m*pow(self.z[4],2)
        self.z[9] = torque*self.z[3] - 0.5*self.z[4]*(force+self.z[8]*self.z[2])
        self.z[12] = self.z[9]/self.z[11]
        omegap = self.z[12]

        # store the results in f and return
        f = zeros_like(x)
        f[0] = omegap
        f[1] = thetap

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
        u[0] = 10*sin(0.5235987755982988+6.283185307179586*T)
        u[1] = 0
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
        # declare the parameters
        i = self.parameters['i']
        m = self.parameters['m']
        l = self.parameters['l']
        g = self.parameters['g']


        # declare the parameters


        # declare the states
        omega = x[0]
        theta = x[1]

        # these are dependent variables that may be needed for the main
        # calculations

        # calculate the outputs
        k = 0.125*(m*pow(self.z[4],2)+4*i*pow(self.z[3],2))*pow(omega,2)
        p = -0.5*g*l*m*self.z[1]
        longoutput = 2 + 2*theta + p + self.z[1] + omega + 0.125*(m*pow(self.z[4],2)+4*i*pow(self.z[3],2))*pow(omega,2)

        # store the results in y and return
        y = zeros(len(self.outputNames))
        y[0] = omega
        y[1] = theta
        y[2] = k
        y[3] = p
        y[4] = longoutput


        return y
