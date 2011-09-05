from numpy import zeros
from numpy import sin, cos, tan
from DynamicSystem import DynamicSystem

class DoublePendulum(DynamicSystem):
    """
<description>

    """

    # model name
    name = 'DoublePendulum'

    filename = ''.join(name.split())
    directory = 'models/' + filename + '/'

    # numerical integration parameters
    intOpts = {'abserr' : 1e-08,
               'relerr' : 1e-07,
               'tf' : 1.0,
               'ti' : 0.0,
               'ts' : 0.1}

    # parameter names and their values
    parameters = {'g' : 9.81,
                  'i1' : 0.5,
                  'i2' : 0.5,
                  'l1' : 2.0,
                  'l2' : 2.0,
                  'm1' : 4.0,
                  'm2' : 4.0}

    # state names
    stateNames = ['omega1',
                  'omega2',
                  'theta1',
                  'theta2']

    # sets the initial conditions of the states
    initialConditions = [0.0,
                         0.0,
                         0.0,
                         0.0]

    # input names
    inputNames = ['torque',
                  'force']

    # output names
    outputNames = ['omega2',
                   'theta2',
                   'kinetic',
                   'potential',
                   'energy']

    # initialize state vector
    x = zeros(len(stateNames))

    # initialize output vector
    y = zeros(len(outputNames))

    # initialize input vector
    u = zeros(len(inputNames))

    # initializes the zees
    z = zeros(67)

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

        self.z[16] = g*m1
        self.z[17] = g*m2
        self.z[63] = l1*m1
        self.z[65] = g*l2*m2

    def f(self, x, t):
        '''Returns the time derivative of the state vector.

        Parameters:
        -----------
        x : ndarray, shape(n,)
            State vector
        t : float
            Time

        Returns:
        --------
        f : ndarray, shape(n,)
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

        # equations of motion
        theta1p = omega1
        theta2p = omega2
        self.z[1] = cos(theta1)
        self.z[2] = sin(theta1)
        self.z[5] = pow(self.z[1],2) + pow(self.z[2],2)
        self.z[3] = cos(theta2)
        self.z[4] = sin(theta2)
        self.z[6] = pow(self.z[3],2) + pow(self.z[4],2)
        self.z[25] = i2*self.z[6]
        self.z[9] = l2*self.z[6]
        self.z[8] = l2*self.z[5]
        self.z[7] = l1*self.z[5]
        self.z[27] = self.z[5]*self.z[25] + 0.25*m2*self.z[9]*(self.z[8]+2*self.z[3]*self.z[7])
        self.z[23] = i1*self.z[5]
        self.z[24] = i2*self.z[5]
        self.z[26] = self.z[5]*self.z[23] + self.z[5]*self.z[24] + 0.25*m1*pow(self.z[7],2) + 0.25*m2*(pow(self.z[8],2)+4*pow(self.z[7],2)+4*self.z[3]*self.z[7]*self.z[8])
        self.z[10] = self.z[5]*omega1
        self.z[11] = self.z[5]*omega1 + self.z[6]*omega2
        self.z[12] = self.z[7]*omega1
        self.z[13] = self.z[10]*self.z[12]
        self.z[14] = -0.5*self.z[8]*omega1 - 0.5*self.z[9]*omega2
        self.z[15] = self.z[11]*self.z[14]
        self.z[19] = self.z[1]*self.z[4] + self.z[2]*self.z[3]
        self.z[28] = m2*self.z[4]*(self.z[8]*self.z[13]+2*self.z[7]*self.z[15])
        self.z[29] = self.z[6]*self.z[24] + 0.25*m2*self.z[9]*(self.z[8]+2*self.z[3]*self.z[7])
        self.z[30] = self.z[6]*self.z[25] + 0.25*m2*pow(self.z[9],2)
        self.z[31] = m2*self.z[4]*self.z[9]*self.z[13]
        self.z[34] = self.z[26]*self.z[30] - self.z[27]*self.z[29]
        self.z[22] = self.z[9]*(force+self.z[17]*self.z[19])
        self.z[33] = 0.5*self.z[22] + 0.5*self.z[31]
        self.z[21] = torque*self.z[5] + force*self.z[4]*self.z[7] - 0.5*force*self.z[8] - self.z[17]*self.z[2]*self.z[7] -force*self.z[3]*self.z[7] - 0.5*self.z[16]*self.z[2]*self.z[7] - 0.5*self.z[17]*self.z[8]*self.z[19]
        self.z[32] = 0.5*self.z[28] - self.z[21]
        self.z[35] = (self.z[27]*self.z[33]-self.z[30]*self.z[32])/self.z[34]
        omega1p = self.z[35]
        self.z[36] = (self.z[26]*self.z[33]-self.z[29]*self.z[32])/self.z[34]
        omega2p = -self.z[36]

        # plug in the derivatives for returning
        f = zeros(len(self.stateNames))
        for i, name in enumerate(self.stateNames):
            exec('f[' + str(i) + '] = ' + name + 'p')

        return f

    def inputs(self, t):
        '''Returns the inputs to the system.

        Parameters:
        -----------
        t : float
            Time

        Returns:
        --------
        u : ndarray, shape(p,)
            Inputs a time t.

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

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
        # defines the parameters locally from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        for i, name in enumerate(self.stateNames):
            exec(name + ' = ' + 'x[' + str(i) + ']')

        # these are dependent variables that may be needed for the main
        # calculations

        # calculate the outputs
        self.z[18] = self.z[1]*self.z[3] - self.z[2]*self.z[4]
        self.z[20] = -self.z[1]*self.z[4] - self.z[2]*self.z[3]
        kinetic = 0.125*m1*pow(self.z[7],2)*pow(omega1,2) + 0.5*i1*pow(self.z[5],2)*pow(omega1,2) + 0.5*i2*pow((self.z[5]*omega1+self.z[6]*omega2),2) + 0.125*m2*(4*pow(self.z[7],2)*pow(omega1,2)+pow((self.z[8]*omega1+self.z[9]*omega2),2)+4*self.z[3]*self.z[7]*omega1*(self.z[8]*omega1+self.z[9]*omega2))
        potential = -0.5*g*(l1*m1*self.z[1]+m2*(l2*self.z[18]+2*l1*self.z[1]))
        energy = potential + kinetic

        # plug in the derivatives for returning
        y = zeros(len(self.outputNames))
        for i, name in enumerate(self.outputNames):
            exec('y[' + str(i) + '] = ' + name)

        return y
