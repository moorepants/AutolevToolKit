from numpy import zeros, zeros_like, dot
from numpy import linspace, rank
from numpy.linalg import eig
from scipy.integrate import odeint
from matplotlib.pyplot import plot, show, legend, xlabel, title, figure
from matplotlib.pyplot import scatter, colorbar, cm, grid, axis
import pickle
import os

class DynamicSystem(object):
    """
    Dynamic System class.

    """

    # model name
    name = 'DynamicSystem'

    filename = ''.join(name.split())
    directory = os.path.join('..', 'models', filename)

    # numerical integration parameters
    intOpts = {'ti' : 0.0,
               'tf' : 1.0,
               'ts' : 0.1,
               'abserr' : 1e-8,
               'relerr' : 1e-07}

    # parameter names and their values
    parameters = {'a' : 1.0,
                  'b' : 2.0}

    # state names
    stateNames = ['x1',
                  'x2']

    # sets the initial conditions of the states
    initialConditions = [0.0,
                         0.0]

    # input names
    inputNames = ['u1']

    # output names
    outputNames = ['y1',
                   'y2']

    # initialize state vector
    x = zeros(len(stateNames))

    # initialize output vector
    y = zeros(len(outputNames))

    # initialize input vector
    u = zeros(len(inputNames))

    # initialize the zees
    z = zeros(1)

    # sets the time to the initial time
    t = intOpts['ti']

    def f(self, x, t):
        '''
        Returns the derivative of the states at the specified time.

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
        a = self.parameters['a']
        b = self.parameters['b']

        # declare the states
        x1 = x[0]
        x2 = x[1]

        # calculates and declare the inputs
        u = self.inputs(t)
        F1 = u[0]

        # calculate the zees
        self.z[0] = 0.

        # calculates the derivatives of the states
        x1p = x2
        x2p = a * 1. + b * F1

        # store the results in f and return
        f = zeros_like(x)
        f[0] = x1p
        f[1] = x2p

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
        # initialize the input array
        u = zeros(len(self.inputNames))
        # calculate or specifiy the input vector
        u[0] = 1.
        return u

    def outputs(self, x):
        '''
        Returns the outputs of the system.

        Parameters
        ----------
        x : ndarray, shape(n,)
            The current state vector.

        Returns
        -------
        y : ndarray, shape(m,)
            The output vector.

        '''
        y = zeros(len(self.outputNames))
        y[0] = x[0]
        y[1] = x[1]

        return y

    def set_initial_conditions(self, *args):
        """Sets a given initial condition.

        This function takes a series of variable, value pairs.

        Parameters
        ----------
        var : string
            The state name of the initial condition.
        val : float
            The value of the initial conidition.

        """

        if len(args) % 2 != 0:
            raise StandardError('You must have an even number of arguments.')

        states = [x for i, x in enumerate(args) if i % 2 == 0]
        values = [x for i, x in enumerate(args) if i % 2 != 0]

        for var, val in zip(states, values):
            try:
                index = self.stateNames.index(var)
            except ValueError:
                raise ValueError('{} is not a valid state.'.format(var))

            self.initialConditions[index] = val

    def set_parameters(self, par):
        """Sets the parameters with the given dictionary.

        Parameters
        ----------
        par : dictionary
            A dictionary at least some or all of the parameters of the model.

        """
        for p in self.parameters.keys():
            try:
                self.parameters[p] = par[p]
            except KeyError:
                print('{}'.format(p) + ' was not in the provided ' +
                        'parameters and thus was not set.')

        # recalculate the constants
        self.constants()

    def simulate(self):
        '''
        Simulates the system.

        Parameters
        ----------

        Returns
        -------

        '''
        # make sure the constants are updated
        self.constants()

        # time vector
        t = linspace(self.intOpts['ti'],
                     self.intOpts['tf'] - self.intOpts['ts'],
                     (self.intOpts['tf'] - self.intOpts['ti']) / self.intOpts['ts'])

        #print self.stateNames
        #print self.outputNames
        #print self.inputNames

        # initialize the vectors
        x = zeros((len(t), len(self.x)))
        u = zeros((len(t), len(self.u)))
        y = zeros((len(t), len(self.y)))

        # set the initial conditions
        x[0] = self.initialConditions
        u[0] = self.inputs(t[0])
        y[0] = self.outputs(x[0])

        for i in range(len(t) - 1):
            print 't =', t[i]
            # set the interval
            t_int = [t[i], t[i + 1]]
            #print "self.t before int = ", self.t
            #print "self.u before int = ", self.u
            #print "self.x before int = ", self.x
            #print "self.z before int = ", self.z
            #print "self.y before int = ", self.y
            # return the next state
            x[i + 1] = odeint(self.f, x[i], t_int)[1, :]
            # calculate the next input value
            u[i + 1] = self.inputs(t[i + 1])
            # calculate the outputs and store them
            y[i + 1] = self.outputs(x[i + 1])
            # update all the attributes
            self.t = t[i + 1]
            self.x = x[i + 1]
            self.u = u[i + 1]
            self.y = y[i + 1]
            #print "self.t after int = ", self.t
            #print "self.u after int = ", self.u
            #print "self.x after int = ", self.x
            #print "self.z after int = ", self.z
            #print "self.y after int = ", self.y

        # make a dictionary of the integration and save it to file
        self.simResults = {'t':t,
                           'x':x,
                           'y':y,
                           'u':u,
                           'model':self.name,
                           'params':self.parameters}

    def save_sim(self):
        '''
        Save simulation to file

        '''
        if os.path.isdir(self.directory):
            pass
        else:
            os.system('mkdir ' + self.directory)
        pickle.dump(self.simResults, open(self.directory + self.filename + '.p', 'w'))

    def plot(self):
        '''
        Makes a plot of the simulation

        '''
        fig = figure()
        plot(self.simResults['t'], self.simResults['y'])
        legend(self.outputNames)
        xlabel('Time [sec]')
        return fig

class LinearDynamicSystem(DynamicSystem):

    name = "LinearDynamicSystem"

    def f(self, x, t):
        '''Returns the derivative of the states'''

        u = self.inputs(t)

        xd = dot(self.A, x) + dot(self.B, u)

        return xd

    def outputs(self, x):

        # the feedforward needs to be included with the inputs, thus making
        # this a function of time.
        y = dot(self.C, x)

        return y

    def linear(self, equilibriumPoint):
        """Calculates the state, input, output and feedforward  matrices for the
        system linearized about the provided equilibrium point.

        Parameters
        ----------
        equilibriumPoint : ndarray, shape(n,)
            The point at which to linearize the system about. The order of the
            values corresponds to the system states.

        """

        self.equilibriumPoint = equilibriumPoint

        # sets the zees for the equilbrium points
        nonlin = DynamicSystem()
        nonlin.f(self.equilibriumPoint, 0.)

        # defines the A, B, C, D matrices
        self.A = zeros((2,2))
        self.B = zeros(2)
        self.C = zeros((5,2))
        self.D = zeros(5)
        self.A[0,0] = 0
        self.A[0,1] = 1
        self.A[1,0] = -0.5*self.z[13]
        self.A[1,1] = 0
        self.B[0] = 0
        self.B[1] = self.z[3]/self.z[11]
        self.C[0,0] = 1
        self.C[0,1] = 0
        self.C[1,0] = 0
        self.C[1,1] = 1
        self.C[2,0] = 0
        self.C[2,1] = 0.25*self.z[14]
        self.C[3,0] = 0.5*self.z[16]
        self.C[3,1] = 0
        self.C[4,0] = 2
        self.C[4,1] = 0
        self.D[0] = 0
        self.D[1] = 0
        self.D[2] = 0
        self.D[3] = 0
        self.D[4] = 0

    def root_loci(self, var, start, stop, num=None):
        """Returns the eigenvalues and eigenvectors as a function of a single
        parameter.

        Parameters
        ----------
        var : string
            The parameter or equilibrim point (state) to vary.
        start : float
            The starting value.
        stop : float
            The ending value.
        num : integer, optional
            The number of steps.

        Returns
        -------
        eValues : ndarray, shape(n,m)
        eVectors : ndarray, shape(n,m,m)
        values : ndarray, shape(n)
            The values at which the model was linea

        """
        eValues = zeros((num, len(self.stateNames)), dtype=complex)
        eVectors = zeros((num, len(self.stateNames), len(self.stateNames)),
                dtype=complex)
        if num is None:
            num = 50
        values = linspace(start, stop, num=num)

        if var in self.parameters.keys():
            original = self.parameters[var]
            for i, val in enumerate(values):
                self.parameters[var] = val
                self.linear(self.equilibriumPoint)
                eValues[i], eVectors[i] = self.eig()
            # set the model back to the default
            self.parameters[var] = original
            self.linear(self.equilibriumPoint)
        elif var in self.stateNames:
            index = self.stateNames.index(var)
            original = self.equilibriumPoint[index]
            for i, val in enumerate(values):
                self.equilibriumPoint[index] = val
                self.linear(self.equilibriumPoint)
                eValues[i], eVectors[i] = self.eig()
            # set the model back to the default
            self.equilibriumPoint[index] = original
            self.linear(self.equilibriumPoint)
        else:
            raise ValueError('{} is not a valid parameter.'.format(var))

        return eValues, eVectors, values

    def eig(self):
        """Returns the eigenvalues and eigenvectors of the system."""
        return eig(self.A)

    def plot_root_loci(self, parameter, start, stop, num=None, axes='complex',
            parts='both'):
        """Returns a plot of the roots with respect to change in a single
        parameter.

        Parameters
        ----------
        parameter : string
            The parameter to vary.
        start : float
            The starting value of the parameter.
        stop : float
            The ending value of the parameter.
        num : integer, optional
            The number of steps in the parameter.
        axes : string
            If `complex` the roots are plotted in the complex plane with a
            colorbar depicting the variation in the parameter. If `parameter`
            the roots are plotting as the real and imaginary parts with respect
            to the parameter variation.
        parts : string, optional
            Can be set to 'both', 'real', or 'imaginary'. This only applies to
            axes='parameter' plots.

        Returns
        -------
        rootLociFig

        """

        # plot a graph with all the outputs
        eValues, eVectors, parValues = self.root_loci(parameter, start, stop, num=num)
        rootLociFig = figure()
        if axes == 'complex':
            x = eValues.real
            y = eValues.imag
            scatter(x, y, s=2, c=parValues, cmap=cm.gist_rainbow, edgecolors='none')
            #colorbar()
            grid()
            axis('equal')
        elif axes == 'parameter':
            if parts == 'both' or parts == 'imaginary':
                plot(parValues, eValues.imag, 'r.')
            if parts == 'both' or parts == 'real':
                plot(parValues, eValues.real, 'k.')
            if start > stop:
                ax = rootLociFig.axes[0]
                ax.invert_xaxis()

        title('Root loci with to {}'.format(parameter))

        return rootLociFig
