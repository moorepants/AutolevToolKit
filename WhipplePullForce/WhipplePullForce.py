from numpy import zeros
from numpy import sin, cos, tan
from DynamicSystem import DynamicSystem

class WhipplePullForce(DynamicSystem):
    """
<description>

    """

    # model name
    name = 'WhipplePullForce'

    filename = ''.join(name.split())
    directory = 'models/' + filename + '/'

    # numerical integration parameters
    intOpts = {'abserr' : 1e-08,
               'relerr' : 1e-07,
               'tf' : 1.0,
               'ti' : 0.0,
               'ts' : 0.1}

    # parameter names and their values
    parameters = {'IBxx' : 9.2,
                  'IBxz' : 2.4,
                  'IByy' : 11.0,
                  'IBzz' : 2.8,
                  'IFxx' : 0.1405,
                  'IFyy' : 0.28,
                  'IHxx' : 0.05892,
                  'IHxz' : -0.00756,
                  'IHyy' : 0.06,
                  'IHzz' : 0.00708,
                  'IRxx' : 0.0603,
                  'IRyy' : 0.12,
                  'c' : 0.08,
                  'g' : 9.81,
                  'lam' : 0.314159265359,
                  'mB' : 85.0,
                  'mF' : 3.0,
                  'mH' : 4.0,
                  'mR' : 2.0,
                  'rF' : 0.35,
                  'rR' : 0.3,
                  'w' : 1.02,
                  'xB' : 0.3,
                  'xH' : 0.9,
                  'xpf' : 0.23,
                  'zB' : -0.9,
                  'zH' : -0.7,
                  'zpf' : -0.91}

    # state names
    stateNames = ['q1',
                  'q2',
                  'q3',
                  'q4',
                  'q5',
                  'q6',
                  'q7',
                  'q8',
                  'u4',
                  'u6',
                  'u7']

    # sets the initial conditions of the states
    initialConditions = [0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         0.0]

    # input names
    inputNames = ['Fphi',
                  'Tdelta',
                  'Tphi',
                  'TthetaR']

    # output names
    outputNames = ['q1',
                   'q2',
                   'q3',
                   'q4',
                   'q5',
                   'q6',
                   'q7',
                   'q8',
                   'u1',
                   'u2',
                   'u3',
                   'u4',
                   'u5',
                   'u6',
                   'u7',
                   'u8',
                   'q9',
                   'q10']

    # initialize state vector
    x = zeros(len(stateNames))

    # initialize output vector
    y = zeros(len(outputNames))

    # initialize input vector
    u = zeros(len(inputNames))

    # initializes the zees
    z = zeros(4232)

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

        self.parameters['d1'] = cos(lam)*(c+w-rR*tan(lam))
        self.parameters['d3'] = -cos(lam)*(c-rF*tan(lam))
        self.parameters['d2'] = -(rF-rR-sin(lam)*self.parameters['d1']-sin(lam)*self.parameters['d3'])/cos(lam)
        self.parameters['id11'] = IRxx
        self.parameters['id22'] = IRyy
        self.parameters['id33'] = IRxx
        self.parameters['ic11'] = IBzz*pow(sin(lam),2) + cos(lam)*(IBxx*cos(lam)-2*IBxz*sin(lam))
        self.parameters['ic12'] = 0
        self.parameters['ic22'] = IByy
        self.parameters['ic23'] = 0
        self.parameters['ic31'] = sin(lam)*(IBxx*cos(lam)-IBxz*sin(lam)) + cos(lam)*(IBxz*cos(lam)-IBzz*sin(lam))
        self.parameters['ic33'] = IBzz*pow(cos(lam),2) + sin(lam)*(IBxx*sin(lam)+2*IBxz*cos(lam))
        self.parameters['ie11'] = IHzz*pow(sin(lam),2) + cos(lam)*(IHxx*cos(lam)-2*IHxz*sin(lam))
        self.parameters['ie12'] = 0
        self.parameters['ie22'] = IHyy
        self.parameters['ie23'] = 0
        self.parameters['ie31'] = sin(lam)*(IHxx*cos(lam)-IHxz*sin(lam)) + cos(lam)*(IHxz*cos(lam)-IHzz*sin(lam))
        self.parameters['ie33'] = IHzz*pow(cos(lam),2) + sin(lam)*(IHxx*sin(lam)+2*IHxz*cos(lam))
        self.parameters['if11'] = IFxx
        self.parameters['if22'] = IFyy
        self.parameters['if33'] = IFxx
        self.parameters['l1'] = xB*cos(lam) - rR*sin(lam) - zB*sin(lam)
        self.parameters['l2'] = rR*cos(lam) + xB*sin(lam) + zB*cos(lam)
        self.parameters['l3'] = xH*cos(lam) - c*cos(lam) - w*cos(lam) - zH*sin(lam)
        self.parameters['l4'] = rR*cos(lam) + xH*sin(lam) + zH*cos(lam)
        self.parameters['massc'] = mB
        self.parameters['massd'] = mR
        self.parameters['masse'] = mH
        self.parameters['massf'] = mF
        self.parameters['l5'] = xpf*cos(lam) - rR*sin(lam) - zpf*sin(lam)
        self.parameters['l6'] = rR*cos(lam) + xpf*sin(lam) + zpf*cos(lam)
        self.z[513] = g*self.parameters['massc']
        self.z[514] = g*self.parameters['massd']
        self.z[515] = g*self.parameters['masse']
        self.z[516] = g*self.parameters['massf']
        self.z[747] = self.parameters['massc'] + self.parameters['massd'] + self.parameters['masse'] + self.parameters['massf']
        self.z[1487] = rR*self.parameters['massd']

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
        self.z[3] = cos(q4)
        self.z[2] = sin(q3)
        self.z[4] = sin(q4)
        self.z[19] = self.z[2]*self.z[4]
        self.z[1] = cos(q3)
        self.z[20] = self.z[1]*self.z[4]
        self.z[38] = pow(self.z[3],2) + pow(self.z[19],2) + pow(self.z[20],2)
        self.z[101] = rR*self.z[38]
        self.z[5] = cos(lam+q5)
        self.z[7] = cos(q6)
        self.z[6] = sin(lam+q5)
        self.z[8] = sin(q6)
        self.z[48] = self.z[5]*self.z[7] - self.z[6]*self.z[8]
        self.z[52] = self.z[38]*self.z[48]
        self.z[49] = self.z[5]*self.z[8] + self.z[6]*self.z[7]
        self.z[54] = self.z[38]*self.z[49]
        self.z[136] = rR*(self.z[48]*self.z[52]+self.z[49]*self.z[54])
        self.z[140] = self.z[101] - self.z[136]
        self.z[208] = self.z[3]*self.z[140]
        self.z[10] = sin(q7)
        self.z[16] = self.z[6]*self.z[10]
        self.z[9] = cos(q7)
        self.z[26] = self.z[3]*self.z[16] + self.z[4]*self.z[9]
        self.z[15] = self.z[5]*self.z[10]
        self.z[70] = self.z[3]*self.z[15]
        self.z[186] = self.z[26]*self.z[70]
        self.z[30] = 1 - pow(self.z[26],2)
        self.z[189] = self.z[186]/pow(self.z[30],0.5)
        self.z[31] = pow(self.z[30],0.5)
        self.z[192] = self.z[189]/pow(self.z[31],2)
        self.z[195] = rF*self.z[192]
        self.z[29] = self.z[3]*self.z[5]
        self.z[77] = self.z[3]*self.z[6]
        self.z[13] = self.z[5]*self.z[9]
        self.z[61] = self.z[3]*self.z[13]
        self.z[198] = (self.z[26]*self.z[189]+self.z[31]*self.z[70])/pow(self.z[31],2)
        self.z[201] = rF*self.z[198]
        self.z[33] = self.z[26]/self.z[31]
        self.z[34] = rF*self.z[33]
        self.z[204] = self.z[195] - d1*self.z[29] - d2*self.z[77] - d3*self.z[61] - self.z[26]*self.z[201] - self.z[34]*self.z[70]
        self.z[17] = self.z[2]*self.z[3]
        self.z[18] = self.z[1]*self.z[3]
        self.z[37] = self.z[17]*self.z[20] - self.z[18]*self.z[19]
        self.z[100] = rR*self.z[37]
        self.z[41] = self.z[1]*self.z[18] + self.z[2]*self.z[17]
        self.z[45] = self.z[5]*self.z[41] + self.z[6]*self.z[37]
        self.z[110] = d1*self.z[45]
        self.z[76] = self.z[4]*self.z[5]
        self.z[40] = self.z[1]*self.z[20] + self.z[2]*self.z[19]
        self.z[113] = d1*self.z[40]
        self.z[14] = self.z[6]*self.z[9]
        self.z[21] = self.z[1]*self.z[13] - self.z[10]*self.z[17] - self.z[14]*self.z[19]
        self.z[73] = self.z[5]*self.z[20] - self.z[2]*self.z[6]
        self.z[22] = self.z[2]*self.z[13] + self.z[10]*self.z[18] + self.z[14]*self.z[20]
        self.z[27] = self.z[1]*self.z[6] + self.z[5]*self.z[19]
        self.z[84] = self.z[21]*self.z[73] + self.z[22]*self.z[27]
        self.z[134] = d3*self.z[84]
        self.z[127] = d2*self.z[84]
        self.z[173] = self.z[4]*self.z[6]
        self.z[174] = self.z[3]*self.z[10] + self.z[9]*self.z[173]
        self.z[25] = self.z[9]*self.z[18] - self.z[2]*self.z[15] - self.z[16]*self.z[20]
        self.z[24] = self.z[16]*self.z[19] - self.z[1]*self.z[15] - self.z[9]*self.z[17]
        self.z[55] = -self.z[2]*self.z[13] - self.z[10]*self.z[18] - self.z[14]*self.z[20]
        self.z[85] = self.z[21]*self.z[25] + self.z[24]*self.z[55]
        self.z[28] = self.z[2]*self.z[6] - self.z[5]*self.z[20]
        self.z[62] = self.z[2]*self.z[15] + self.z[16]*self.z[20] - self.z[9]*self.z[18]
        self.z[78] = self.z[24]*self.z[28] + self.z[27]*self.z[62]
        self.z[128] = d3*self.z[85] - d2*self.z[78]
        self.z[176] = self.z[3]*self.z[9] - self.z[10]*self.z[173]
        self.z[11] = cos(q8)
        self.z[12] = sin(q8)
        self.z[94] = self.z[11]*self.z[85] + self.z[12]*self.z[78]
        self.z[161] = self.z[34]*self.z[94]
        self.z[179] = self.z[11]*self.z[174] + self.z[12]*self.z[76]
        self.z[32] = 1/self.z[31]
        self.z[35] = rF*self.z[32]
        self.z[90] = self.z[11]*self.z[78] - self.z[12]*self.z[85]
        self.z[141] = self.z[11]*self.z[21] - self.z[12]*self.z[27]
        self.z[144] = self.z[11]*self.z[27] + self.z[12]*self.z[21]
        self.z[153] = self.z[35]*(self.z[24]*self.z[84]+self.z[90]*self.z[141]+self.z[94]*self.z[144])
        self.z[142] = self.z[11]*self.z[22] - self.z[12]*self.z[28]
        self.z[145] = self.z[11]*self.z[28] + self.z[12]*self.z[22]
        self.z[158] = self.z[35]*(self.z[25]*self.z[84]+self.z[90]*self.z[142]+self.z[94]*self.z[145])
        self.z[147] = self.z[34]*self.z[90]
        self.z[182] = self.z[12]*self.z[174] - self.z[11]*self.z[76]
        self.z[214] = self.z[3]*self.z[100] + self.z[3]*self.z[110] + self.z[76]*self.z[113] + self.z[76]*self.z[134] + self.z[127]*self.z[174] + self.z[128]*self.z[176] + self.z[161]*self.z[179] - self.z[1]*self.z[153] - self.z[2]*self.z[158] - self.z[147]*self.z[182]
        self.z[239] = self.z[204]*self.z[214]
        self.z[89] = pow(self.z[11],2) + pow(self.z[12],2)
        self.z[159] = self.z[25]*self.z[35]*self.z[89]
        self.z[154] = self.z[24]*self.z[35]*self.z[89]
        self.z[213] = self.z[1]*self.z[159] - self.z[2]*self.z[154]
        self.z[254] = self.z[2]*self.z[213]
        self.z[218] = -self.z[1]*self.z[154] - self.z[2]*self.z[159]
        self.z[236] = self.z[204]*self.z[218]
        self.z[178] = self.z[11]*self.z[13] - self.z[6]*self.z[12]
        self.z[99] = rR*self.z[40]
        self.z[181] = self.z[6]*self.z[11] + self.z[12]*self.z[13]
        self.z[211] = self.z[1]*self.z[158] + self.z[13]*self.z[127] + self.z[161]*self.z[178] - self.z[99] - self.z[2]*self.z[153] -self.z[6]*self.z[113] - self.z[6]*self.z[134] - self.z[15]*self.z[128] - self.z[147]*self.z[181]
        self.z[252] = self.z[2]*self.z[211]
        self.z[257] = self.z[239]*self.z[254] - self.z[236]*self.z[252]
        self.z[72] = self.z[5]*self.z[17]
        self.z[75] = self.z[5]*self.z[18]
        self.z[23] = self.z[4]*self.z[10] - self.z[3]*self.z[14]
        self.z[83] = self.z[21]*self.z[72] - self.z[22]*self.z[75] - self.z[23]*self.z[76]
        self.z[64] = self.z[9]*self.z[19] + self.z[16]*self.z[17]
        self.z[66] = -self.z[9]*self.z[20] - self.z[16]*self.z[18]
        self.z[68] = self.z[3]*self.z[9] - self.z[4]*self.z[16]
        self.z[80] = self.z[27]*self.z[64] + self.z[28]*self.z[66] + self.z[29]*self.z[68]
        self.z[57] = self.z[10]*self.z[19] - self.z[14]*self.z[17]
        self.z[58] = self.z[14]*self.z[18] - self.z[10]*self.z[20]
        self.z[60] = self.z[3]*self.z[10] + self.z[4]*self.z[14]
        self.z[87] = self.z[24]*self.z[57] + self.z[25]*self.z[58] + self.z[26]*self.z[60]
        self.z[92] = self.z[11]*self.z[80] - self.z[12]*self.z[87]
        self.z[96] = self.z[11]*self.z[87] + self.z[12]*self.z[80]
        self.z[157] = self.z[35]*(self.z[25]*self.z[83]+self.z[92]*self.z[142]+self.z[96]*self.z[145])
        self.z[126] = d2*self.z[83]
        self.z[163] = self.z[34]*self.z[96]
        self.z[39] = self.z[1]*self.z[17] - self.z[2]*self.z[18]
        self.z[98] = rR*self.z[39]
        self.z[152] = self.z[35]*(self.z[24]*self.z[83]+self.z[92]*self.z[141]+self.z[96]*self.z[144])
        self.z[112] = d1*self.z[39]
        self.z[133] = d3*self.z[83]
        self.z[130] = d3*self.z[87] - d2*self.z[80]
        self.z[149] = self.z[34]*self.z[92]
        self.z[210] = self.z[1]*self.z[157] + self.z[13]*self.z[126] + self.z[163]*self.z[178] - self.z[98] - self.z[2]*self.z[152] -self.z[6]*self.z[112] - self.z[6]*self.z[133] - self.z[15]*self.z[130] - self.z[149]*self.z[181]
        self.z[50] = -self.z[5]*self.z[8] - self.z[6]*self.z[7]
        self.z[51] = self.z[37]*self.z[48] + self.z[41]*self.z[50]
        self.z[53] = self.z[37]*self.z[49] + self.z[41]*self.z[48]
        self.z[135] = rR*(self.z[48]*self.z[51]+self.z[49]*self.z[53])
        self.z[139] = self.z[100] - self.z[135]
        self.z[207] = self.z[3]*self.z[139]
        self.z[228] = self.z[2]*self.z[207]
        self.z[184] = self.z[26]*self.z[68]
        self.z[187] = self.z[184]/pow(self.z[30],0.5)
        self.z[190] = self.z[187]/pow(self.z[31],2)
        self.z[193] = rF*self.z[190]
        self.z[196] = (self.z[26]*self.z[187]+self.z[31]*self.z[68])/pow(self.z[31],2)
        self.z[199] = rF*self.z[196]
        self.z[202] = self.z[193] + rR*self.z[4] + d1*self.z[173] + d3*self.z[60] - d2*self.z[76] - self.z[26]*self.z[199] -self.z[34]*self.z[68]
        self.z[221] = self.z[211]*self.z[218] - self.z[213]*self.z[214]
        self.z[42] = pow(self.z[5],2) + pow(self.z[6],2)
        self.z[137] = rR*self.z[42]
        self.z[206] = self.z[3]*self.z[18] + self.z[4]*self.z[20]
        self.z[231] = self.z[137]*self.z[206]
        self.z[71] = self.z[1]*self.z[5] - self.z[6]*self.z[19]
        self.z[74] = self.z[2]*self.z[5] + self.z[6]*self.z[20]
        self.z[82] = self.z[21]*self.z[71] + self.z[22]*self.z[74] - self.z[23]*self.z[77]
        self.z[65] = self.z[1]*self.z[16] + self.z[15]*self.z[19]
        self.z[67] = self.z[2]*self.z[16] - self.z[15]*self.z[20]
        self.z[81] = self.z[27]*self.z[65] + self.z[28]*self.z[67] + self.z[29]*self.z[70]
        self.z[56] = -self.z[1]*self.z[14] - self.z[13]*self.z[19]
        self.z[59] = self.z[13]*self.z[20] - self.z[2]*self.z[14]
        self.z[86] = self.z[24]*self.z[56] + self.z[25]*self.z[59] - self.z[26]*self.z[61]
        self.z[93] = self.z[11]*self.z[81] - self.z[12]*self.z[86]
        self.z[95] = self.z[11]*self.z[86] + self.z[12]*self.z[81]
        self.z[156] = self.z[35]*(self.z[25]*self.z[82]+self.z[93]*self.z[142]+self.z[95]*self.z[145])
        self.z[125] = d2*self.z[82]
        self.z[162] = self.z[34]*self.z[95]
        self.z[151] = self.z[35]*(self.z[24]*self.z[82]+self.z[93]*self.z[141]+self.z[95]*self.z[144])
        self.z[114] = d1*self.z[42]
        self.z[132] = d3*self.z[82]
        self.z[129] = d3*self.z[86] - d2*self.z[81]
        self.z[150] = self.z[34]*self.z[93]
        self.z[209] = self.z[1]*self.z[156] + self.z[13]*self.z[125] + self.z[162]*self.z[178] - self.z[2]*self.z[151] - self.z[6]*self.z[114] - self.z[6]*self.z[132] - self.z[15]*self.z[129] - self.z[150]*self.z[181]
        self.z[216] = self.z[76]*self.z[114] + self.z[76]*self.z[132] + self.z[125]*self.z[174] + self.z[129]*self.z[176] +self.z[162]*self.z[179] - self.z[1]*self.z[151] - self.z[2]*self.z[156] - self.z[150]*self.z[182]
        self.z[222] = self.z[209]*self.z[218] - self.z[213]*self.z[216]
        self.z[225] = self.z[2]*self.z[218] - self.z[1]*self.z[213]
        self.z[233] = self.z[137]*self.z[207]
        self.z[268] = self.z[221]*self.z[231] + self.z[222]*self.z[228] - self.z[225]*self.z[233]
        self.z[46] = self.z[6]*self.z[38]
        self.z[111] = d1*self.z[46]
        self.z[215] = self.z[3]*self.z[101] + self.z[3]*self.z[111] + self.z[76]*self.z[112] + self.z[76]*self.z[133] + self.z[126]*self.z[174] + self.z[130]*self.z[176] + self.z[163]*self.z[179] - self.z[1]*self.z[152] - self.z[2]*self.z[157] - self.z[149]*self.z[182]
        self.z[265] = self.z[204]*self.z[213]
        self.z[223] = self.z[1]*self.z[207]
        self.z[205] = -self.z[3]*self.z[17] - self.z[4]*self.z[19]
        self.z[219] = self.z[1]*self.z[206] - self.z[2]*self.z[205]
        self.z[230] = self.z[1]*self.z[218] + self.z[2]*self.z[213]
        self.z[234] = self.z[204]*(self.z[223]*self.z[225]-self.z[219]*self.z[221]-self.z[228]*self.z[230])
        self.z[274] = (self.z[208]*self.z[257]+self.z[210]*self.z[228]*self.z[236]-self.z[202]*self.z[268]-self.z[215]*self.z[228]*self.z[265])/self.z[234]
        self.z[63] = self.z[10]*self.z[17] + self.z[14]*self.z[19] - self.z[1]*self.z[13]
        self.z[69] = self.z[3]*self.z[14] - self.z[4]*self.z[10]
        self.z[79] = self.z[27]*self.z[63] + self.z[28]*self.z[55] + self.z[29]*self.z[69]
        self.z[88] = pow(self.z[24],2) + pow(self.z[25],2) + pow(self.z[26],2)
        self.z[91] = self.z[11]*self.z[79] - self.z[12]*self.z[88]
        self.z[97] = self.z[11]*self.z[88] + self.z[12]*self.z[79]
        self.z[160] = self.z[35]*(self.z[91]*self.z[142]+self.z[97]*self.z[145])
        self.z[164] = self.z[34]*self.z[97]
        self.z[155] = self.z[35]*(self.z[91]*self.z[141]+self.z[97]*self.z[144])
        self.z[131] = d3*self.z[88] - d2*self.z[79]
        self.z[148] = self.z[34]*self.z[91]
        self.z[212] = self.z[1]*self.z[160] + self.z[164]*self.z[178] - self.z[2]*self.z[155] - self.z[15]*self.z[131] - self.z[148]*self.z[181]
        self.z[185] = self.z[26]*self.z[69]
        self.z[188] = self.z[185]/pow(self.z[30],0.5)
        self.z[191] = self.z[188]/pow(self.z[31],2)
        self.z[194] = rF*self.z[191]
        self.z[197] = (self.z[26]*self.z[188]+self.z[31]*self.z[69])/pow(self.z[31],2)
        self.z[200] = rF*self.z[197]
        self.z[203] = self.z[194] + d3*self.z[26] - self.z[26]*self.z[200] - self.z[34]*self.z[69]
        self.z[217] = self.z[131]*self.z[176] + self.z[164]*self.z[179] - self.z[1]*self.z[155] - self.z[2]*self.z[160] - self.z[148]*self.z[182]
        self.z[275] = (self.z[212]*self.z[228]*self.z[236]-self.z[203]*self.z[268]-self.z[217]*self.z[228]*self.z[265])/self.z[234]
        self.z[47] = pow(self.z[7],2) + pow(self.z[8],2)
        self.z[138] = rR*self.z[47]
        self.z[238] = self.z[206]*self.z[213]
        self.z[235] = self.z[206]*self.z[211] - self.z[2]*self.z[207]
        self.z[241] = self.z[207]*self.z[213]
        self.z[242] = self.z[1]*self.z[204]
        self.z[243] = self.z[238]*self.z[239] - self.z[235]*self.z[236] - self.z[241]*self.z[242]
        self.z[273] = self.z[138]*self.z[243]/self.z[234]
        u1 = self.z[274]*u4 + self.z[275]*u7 - self.z[273]*u6
        q1p = u1
        self.z[246] = self.z[205]*self.z[213]
        self.z[247] = self.z[2]*self.z[204]
        self.z[244] = self.z[205]*self.z[211] - self.z[1]*self.z[207]
        self.z[248] = self.z[239]*self.z[246] + self.z[241]*self.z[247] - self.z[236]*self.z[244]
        self.z[276] = self.z[138]*self.z[248]/self.z[234]
        self.z[260] = self.z[1]*self.z[213]
        self.z[258] = self.z[1]*self.z[211]
        self.z[261] = self.z[239]*self.z[260] - self.z[236]*self.z[258]
        self.z[226] = self.z[137]*self.z[205]
        self.z[269] = self.z[221]*self.z[226] + self.z[222]*self.z[223] - self.z[230]*self.z[233]
        self.z[277] = (self.z[208]*self.z[261]+self.z[210]*self.z[223]*self.z[236]-self.z[202]*self.z[269]-self.z[215]*self.z[223]*self.z[265])/self.z[234]
        self.z[278] = (self.z[212]*self.z[223]*self.z[236]-self.z[203]*self.z[269]-self.z[217]*self.z[223]*self.z[265])/self.z[234]
        u2 = self.z[276]*u6 - self.z[277]*u4 - self.z[278]*u7
        q2p = u2
        self.z[262] = self.z[242]*self.z[260] + self.z[247]*self.z[254]
        self.z[270] = self.z[219]*self.z[222] + self.z[225]*self.z[226] - self.z[230]*self.z[231]
        self.z[280] = (self.z[208]*self.z[262]+self.z[210]*self.z[219]*self.z[236]-self.z[202]*self.z[270]-self.z[215]*self.z[219]*self.z[265])/self.z[234]
        self.z[281] = (self.z[212]*self.z[219]*self.z[236]-self.z[203]*self.z[270]-self.z[217]*self.z[219]*self.z[265])/self.z[234]
        self.z[249] = self.z[2]*self.z[205] - self.z[1]*self.z[206]
        self.z[250] = self.z[238]*self.z[247] + self.z[242]*self.z[246] - self.z[236]*self.z[249]
        self.z[279] = self.z[138]*self.z[250]/self.z[234]
        u3 = self.z[280]*u4 + self.z[281]*u7 - self.z[279]*u6
        q3p = u3
        q4p = u4
        self.z[271] = self.z[219]*self.z[221] + self.z[228]*self.z[230] - self.z[223]*self.z[225]
        self.z[282] = self.z[202]*self.z[271]/self.z[234]
        self.z[283] = self.z[203]*self.z[271]/self.z[234]
        u5 = self.z[282]*u4 + self.z[283]*u7
        q5p = u5
        q6p = u6
        q7p = u7
        self.z[251] = self.z[239]*self.z[249] - self.z[235]*self.z[247] - self.z[242]*self.z[244]
        self.z[284] = self.z[138]*self.z[251]/self.z[234]
        self.z[220] = self.z[211]*self.z[216] - self.z[209]*self.z[214]
        self.z[229] = self.z[1]*self.z[216] + self.z[2]*self.z[209]
        self.z[232] = self.z[1]*self.z[214] + self.z[2]*self.z[211]
        self.z[36] = pow(self.z[1],2) + pow(self.z[2],2)
        self.z[224] = self.z[2]*self.z[216] - self.z[1]*self.z[209]
        self.z[227] = self.z[2]*self.z[214] - self.z[1]*self.z[211]
        self.z[272] = self.z[219]*self.z[220] + self.z[228]*self.z[229] + self.z[231]*self.z[232] - self.z[36]*self.z[233] -self.z[223]*self.z[224] - self.z[226]*self.z[227]
        self.z[264] = self.z[219]*self.z[239] - self.z[223]*self.z[242] - self.z[228]*self.z[247]
        self.z[263] = -self.z[242]*self.z[258] - self.z[247]*self.z[252]
        self.z[266] = self.z[204]*self.z[211]
        self.z[267] = self.z[219]*self.z[266] + self.z[228]*self.z[242] - self.z[223]*self.z[247]
        self.z[285] = (self.z[202]*self.z[272]+self.z[210]*self.z[264]-self.z[208]*self.z[263]-self.z[215]*self.z[267])/self.z[234]
        self.z[286] = (self.z[203]*self.z[272]+self.z[212]*self.z[264]-self.z[217]*self.z[267])/self.z[234]
        u8 = -self.z[284]*u6 - self.z[285]*u4 - self.z[286]*u7
        q8p = u8
        self.z[649] = id22*self.z[47]
        self.z[729] = self.z[39]*self.z[649]
        self.z[757] = self.z[40]*self.z[649]
        self.z[764] = self.z[42]*self.z[649]
        self.z[776] = self.z[729] + self.z[280]*self.z[757] + self.z[282]*self.z[764]
        self.z[116] = l4*self.z[83]
        self.z[120] = l3*self.z[87] - l4*self.z[80]
        self.z[123] = l3*self.z[83]
        self.z[107] = l1*self.z[39]
        self.z[44] = self.z[5]*self.z[38]
        self.z[106] = l1*self.z[46] - l2*self.z[44]
        self.z[102] = l2*self.z[39]
        self.z[732] = masse*(self.z[18]*self.z[101]+self.z[18]*self.z[111]+self.z[22]*self.z[116]+self.z[25]*self.z[120]-self.z[2]*self.z[98]-self.z[28]*self.z[112]-self.z[28]*self.z[123]) + massf*(self.z[18]*self.z[101]+self.z[18]*self.z[111]+self.z[22]*self.z[126]+self.z[25]*self.z[130]-self.z[2]*self.z[98]-self.z[28]*self.z[112]-self.z[28]*self.z[133]) - massd*(self.z[2]*self.z[98]-self.z[18]*self.z[101]) - massc*(self.z[2]*self.z[98]+self.z[28]*self.z[107]-self.z[18]*self.z[101]-self.z[18]*self.z[106]-self.z[74]*self.z[102])
        self.z[117] = l4*self.z[84]
        self.z[118] = l3*self.z[85] - l4*self.z[78]
        self.z[124] = l3*self.z[84]
        self.z[108] = l1*self.z[40]
        self.z[43] = self.z[5]*self.z[37] - self.z[6]*self.z[41]
        self.z[105] = l1*self.z[45] - l2*self.z[43]
        self.z[103] = l2*self.z[40]
        self.z[751] = masse*(self.z[18]*self.z[100]+self.z[18]*self.z[110]+self.z[22]*self.z[117]+self.z[25]*self.z[118]-self.z[2]*self.z[99]-self.z[28]*self.z[113]-self.z[28]*self.z[124]) + massf*(self.z[18]*self.z[100]+self.z[18]*self.z[110]+self.z[22]*self.z[127]+self.z[25]*self.z[128]-self.z[2]*self.z[99]-self.z[28]*self.z[113]-self.z[28]*self.z[134]) - massd*(self.z[2]*self.z[99]-self.z[18]*self.z[100]) - massc*(self.z[2]*self.z[99]+self.z[28]*self.z[108]-self.z[18]*self.z[100]-self.z[18]*self.z[105]-self.z[74]*self.z[103])
        self.z[115] = l4*self.z[82]
        self.z[119] = l3*self.z[86] - l4*self.z[81]
        self.z[122] = l3*self.z[82]
        self.z[109] = l1*self.z[42]
        self.z[104] = l2*self.z[42]
        self.z[752] = masse*(self.z[22]*self.z[115]+self.z[25]*self.z[119]-self.z[28]*self.z[114]-self.z[28]*self.z[122]) +massf*(self.z[22]*self.z[125]+self.z[25]*self.z[129]-self.z[28]*self.z[114]-self.z[28]*self.z[132]) - massc*(self.z[28]*self.z[109]-self.z[74]*self.z[104])
        self.z[779] = self.z[732] + self.z[280]*self.z[751] + self.z[282]*self.z[752] - self.z[747]*self.z[277]
        self.z[733] = masse*(self.z[21]*self.z[116]+self.z[24]*self.z[120]-self.z[1]*self.z[98]-self.z[17]*self.z[101]-self.z[17]*self.z[111]-self.z[27]*self.z[112]-self.z[27]*self.z[123]) + massf*(self.z[21]*self.z[126]+self.z[24]*self.z[130]-self.z[1]*self.z[98]-self.z[17]*self.z[101]-self.z[17]*self.z[111]-self.z[27]*self.z[112]-self.z[27]*self.z[133]) - massd*(self.z[1]*self.z[98]+self.z[17]*self.z[101]) - massc*(self.z[1]*self.z[98]+self.z[17]*self.z[101]+self.z[17]*self.z[106]+self.z[27]*self.z[107]-self.z[71]*self.z[102])
        self.z[748] = masse*(self.z[21]*self.z[117]+self.z[24]*self.z[118]-self.z[1]*self.z[99]-self.z[17]*self.z[100]-self.z[17]*self.z[110]-self.z[27]*self.z[113]-self.z[27]*self.z[124]) + massf*(self.z[21]*self.z[127]+self.z[24]*self.z[128]-self.z[1]*self.z[99]-self.z[17]*self.z[100]-self.z[17]*self.z[110]-self.z[27]*self.z[113]-self.z[27]*self.z[134]) - massd*(self.z[1]*self.z[99]+self.z[17]*self.z[100]) - massc*(self.z[1]*self.z[99]+self.z[17]*self.z[100]+self.z[17]*self.z[105]+self.z[27]*self.z[108]-self.z[71]*self.z[103])
        self.z[749] = masse*(self.z[21]*self.z[115]+self.z[24]*self.z[119]-self.z[27]*self.z[114]-self.z[27]*self.z[122]) +massf*(self.z[21]*self.z[125]+self.z[24]*self.z[129]-self.z[27]*self.z[114]-self.z[27]*self.z[132]) - massc*(self.z[27]*self.z[109]-self.z[71]*self.z[104])
        self.z[780] = self.z[733] + self.z[747]*self.z[274] + self.z[280]*self.z[748] + self.z[282]*self.z[749]
        self.z[610] = ic12*self.z[43]
        self.z[621] = ic22*self.z[40]
        self.z[631] = ic23*self.z[45]
        self.z[647] = id22*self.z[40]
        self.z[607] = ic11*self.z[43]
        self.z[617] = ic12*self.z[40]
        self.z[628] = ic31*self.z[45]
        self.z[613] = ic31*self.z[43]
        self.z[625] = ic23*self.z[40]
        self.z[634] = ic33*self.z[45]
        self.z[643] = id11*self.z[51]
        self.z[651] = id33*self.z[53]
        self.z[660] = ie11*self.z[78]
        self.z[677] = ie12*self.z[84]
        self.z[687] = ie31*self.z[85]
        self.z[665] = ie12*self.z[78]
        self.z[681] = ie22*self.z[84]
        self.z[692] = ie23*self.z[85]
        self.z[715] = if22*self.z[84]
        self.z[670] = ie31*self.z[78]
        self.z[685] = ie23*self.z[84]
        self.z[697] = ie33*self.z[85]
        self.z[708] = if11*self.z[90]
        self.z[718] = if33*self.z[94]
        self.z[726] = self.z[39]*self.z[610] + self.z[39]*self.z[621] + self.z[39]*self.z[631] + self.z[39]*self.z[647] + self.z[44]*self.z[607] + self.z[44]*self.z[617] + self.z[44]*self.z[628] + self.z[46]*self.z[613] + self.z[46]*self.z[625] + self.z[46]*self.z[634] + self.z[52]*self.z[643] + self.z[54]*self.z[651] + self.z[80]*self.z[660] + self.z[80]*self.z[677] + self.z[80]*self.z[687] + self.z[83]*self.z[665] + self.z[83]*self.z[681] + self.z[83]*self.z[692] + self.z[83]*self.z[715] + self.z[87]*self.z[670] + self.z[87]*self.z[685] + self.z[87]*self.z[697] + self.z[92]*self.z[708] + self.z[96]*self.z[718] + massd*(self.z[98]*self.z[99]+self.z[100]*self.z[101]) - massc*(self.z[5]*self.z[98]*self.z[103]+self.z[5]*self.z[99]*self.z[102]-self.z[98]*self.z[99]-self.z[100]*self.z[101]-self.z[100]*self.z[106]-self.z[101]*self.z[105]-self.z[102]*self.z[103]-self.z[105]*self.z[106]-self.z[107]*self.z[108]-self.z[6]*self.z[98]*self.z[108]-self.z[6]*self.z[99]*self.z[107]) - masse*(self.z[13]*self.z[98]*self.z[117]+self.z[13]*self.z[99]*self.z[116]-self.z[98]*self.z[99]-self.z[100]*self.z[101]-self.z[100]*self.z[111]-self.z[101]*self.z[110]-self.z[110]*self.z[111]-self.z[112]*self.z[113]-self.z[112]*self.z[124]-self.z[113]*self.z[123]-self.z[116]*self.z[117]-self.z[118]*self.z[120]-self.z[123]*self.z[124]-self.z[6]*self.z[98]*self.z[113]-self.z[6]*self.z[98]*self.z[124]-self.z[6]*self.z[99]*self.z[112]-self.z[6]*self.z[99]*self.z[123]-self.z[9]*self.z[100]*self.z[120]-self.z[9]*self.z[101]*self.z[118]-self.z[9]*self.z[110]*self.z[120]-self.z[9]*self.z[111]*self.z[118]-self.z[10]*self.z[100]*self.z[116]-self.z[10]*self.z[101]*self.z[117]-self.z[10]*self.z[110]*self.z[116]-self.z[10]*self.z[111]*self.z[117]-self.z[15]*self.z[98]*self.z[118]-self.z[15]*self.z[99]*self.z[120]) - massf*(self.z[13]*self.z[98]*self.z[127]+self.z[13]*self.z[99]*self.z[126]-self.z[98]*self.z[99]-self.z[100]*self.z[101]-self.z[100]*self.z[111]-self.z[101]*self.z[110]-self.z[110]*self.z[111]-self.z[112]*self.z[113]-self.z[112]*self.z[134]-self.z[113]*self.z[133]-self.z[126]*self.z[127]-self.z[128]*self.z[130]-self.z[133]*self.z[134]-self.z[6]*self.z[98]*self.z[113]-self.z[6]*self.z[98]*self.z[134]-self.z[6]*self.z[99]*self.z[112]-self.z[6]*self.z[99]*self.z[133]-self.z[9]*self.z[100]*self.z[130]-self.z[9]*self.z[101]*self.z[128]-self.z[9]*self.z[110]*self.z[130]-self.z[9]*self.z[111]*self.z[128]-self.z[10]*self.z[100]*self.z[126]-self.z[10]*self.z[101]*self.z[127]-self.z[10]*self.z[110]*self.z[126]-self.z[10]*self.z[111]*self.z[127]-self.z[15]*self.z[98]*self.z[128]-self.z[15]*self.z[99]*self.z[130])
        self.z[754] = self.z[40]*self.z[610] + self.z[40]*self.z[621] + self.z[40]*self.z[631] + self.z[40]*self.z[647] + self.z[43]*self.z[607] + self.z[43]*self.z[617] + self.z[43]*self.z[628] + self.z[45]*self.z[613] + self.z[45]*self.z[625] + self.z[45]*self.z[634] + self.z[51]*self.z[643] + self.z[53]*self.z[651] + self.z[78]*self.z[660] + self.z[78]*self.z[677] + self.z[78]*self.z[687] + self.z[84]*self.z[665] + self.z[84]*self.z[681] + self.z[84]*self.z[692] + self.z[84]*self.z[715] + self.z[85]*self.z[670] + self.z[85]*self.z[685] + self.z[85]*self.z[697] + self.z[90]*self.z[708] + self.z[94]*self.z[718] + massd*(pow(self.z[99],2)+pow(self.z[100],2)) - massc*(2*self.z[5]*self.z[99]*self.z[103]-2*self.z[100]*self.z[105]-pow(self.z[99],2)-pow(self.z[100],2)-pow(self.z[103],2)-pow(self.z[105],2)-pow(self.z[108],2)-2*self.z[6]*self.z[99]*self.z[108]) - masse*(2*self.z[13]*self.z[99]*self.z[117]-2*self.z[100]*self.z[110]-2*self.z[113]*self.z[124]-pow(self.z[99],2)-pow(self.z[100],2)-pow(self.z[110],2)-pow(self.z[113],2)-pow(self.z[117],2)-pow(self.z[118],2)-pow(self.z[124],2)-2*self.z[6]*self.z[99]*self.z[113]-2*self.z[6]*self.z[99]*self.z[124]-2*self.z[9]*self.z[100]*self.z[118]-2*self.z[9]*self.z[110]*self.z[118]-2*self.z[10]*self.z[100]*self.z[117]-2*self.z[10]*self.z[110]*self.z[117]-2*self.z[15]*self.z[99]*self.z[118]) - massf*(2*self.z[13]*self.z[99]*self.z[127]-2*self.z[100]*self.z[110]-2*self.z[113]*self.z[134]-pow(self.z[99],2)-pow(self.z[100],2)-pow(self.z[110],2)-pow(self.z[113],2)-pow(self.z[127],2)-pow(self.z[128],2)-pow(self.z[134],2)-2*self.z[6]*self.z[99]*self.z[113]-2*self.z[6]*self.z[99]*self.z[134]-2*self.z[9]*self.z[100]*self.z[128]-2*self.z[9]*self.z[110]*self.z[128]-2*self.z[10]*self.z[100]*self.z[127]-2*self.z[10]*self.z[110]*self.z[127]-2*self.z[15]*self.z[99]*self.z[128])
        self.z[761] = self.z[42]*self.z[647] + self.z[81]*self.z[660] + self.z[81]*self.z[677] + self.z[81]*self.z[687] + self.z[82]*self.z[665] + self.z[82]*self.z[681] + self.z[82]*self.z[692] + self.z[82]*self.z[715] + self.z[86]*self.z[670] + self.z[86]*self.z[685] + self.z[86]*self.z[697] + self.z[93]*self.z[708] + self.z[95]*self.z[718] + self.z[42]*(self.z[610]+self.z[621]+self.z[631]) - massc*(self.z[5]*self.z[99]*self.z[104]-self.z[103]*self.z[104]-self.z[108]*self.z[109]-self.z[6]*self.z[99]*self.z[109]) - masse*(self.z[13]*self.z[99]*self.z[115]-self.z[113]*self.z[114]-self.z[113]*self.z[122]-self.z[114]*self.z[124]-self.z[115]*self.z[117]-self.z[118]*self.z[119]-self.z[122]*self.z[124]-self.z[6]*self.z[99]*self.z[114]-self.z[6]*self.z[99]*self.z[122]-self.z[9]*self.z[100]*self.z[119]-self.z[9]*self.z[110]*self.z[119]-self.z[10]*self.z[100]*self.z[115]-self.z[10]*self.z[110]*self.z[115]-self.z[15]*self.z[99]*self.z[119]) - massf*(self.z[13]*self.z[99]*self.z[125]-self.z[113]*self.z[114]-self.z[113]*self.z[132]-self.z[114]*self.z[134]-self.z[125]*self.z[127]-self.z[128]*self.z[129]-self.z[132]*self.z[134]-self.z[6]*self.z[99]*self.z[114]-self.z[6]*self.z[99]*self.z[132]-self.z[9]*self.z[100]*self.z[129]-self.z[9]*self.z[110]*self.z[129]-self.z[10]*self.z[100]*self.z[125]-self.z[10]*self.z[110]*self.z[125]-self.z[15]*self.z[99]*self.z[129])
        self.z[770] = self.z[89]*self.z[715]
        self.z[773] = self.z[726] + self.z[274]*self.z[748] + self.z[280]*self.z[754] + self.z[282]*self.z[761] - self.z[277]*self.z[751] - self.z[285]*self.z[770]
        self.z[716] = if22*self.z[89]
        self.z[731] = self.z[83]*self.z[716]
        self.z[759] = self.z[84]*self.z[716]
        self.z[766] = self.z[82]*self.z[716]
        self.z[771] = self.z[89]*self.z[716]
        self.z[778] = self.z[731] + self.z[280]*self.z[759] + self.z[282]*self.z[766] - self.z[285]*self.z[771]
        self.z[783] = self.z[776] + self.z[276]*self.z[779] - self.z[273]*self.z[780] - self.z[279]*self.z[773] - self.z[284]*self.z[778]
        self.z[121] = l3*self.z[88] - l4*self.z[79]
        self.z[745] = self.z[25]*(masse*self.z[121]+massf*self.z[131])
        self.z[744] = self.z[24]*(masse*self.z[121]+massf*self.z[131])
        self.z[661] = ie11*self.z[79]
        self.z[690] = ie31*self.z[88]
        self.z[666] = ie12*self.z[79]
        self.z[695] = ie23*self.z[88]
        self.z[671] = ie31*self.z[79]
        self.z[700] = ie33*self.z[88]
        self.z[709] = if11*self.z[91]
        self.z[721] = if33*self.z[97]
        self.z[758] = self.z[78]*self.z[661] + self.z[78]*self.z[690] + self.z[84]*self.z[666] + self.z[84]*self.z[695] + self.z[85]*self.z[671] + self.z[85]*self.z[700] + self.z[90]*self.z[709] + self.z[94]*self.z[721] + masse*self.z[121]*(self.z[118]+self.z[9]*self.z[100]+self.z[9]*self.z[110]+self.z[15]*self.z[99]) + massf*self.z[131]*(self.z[128]+self.z[9]*self.z[100]+self.z[9]*self.z[110]+self.z[15]*self.z[99])
        self.z[791] = self.z[276]*self.z[745] - self.z[273]*self.z[744] - self.z[279]*self.z[758]
        self.z[790] = -self.z[747]*self.z[273] - self.z[279]*self.z[748]
        self.z[736] = self.z[47]*self.z[647]
        self.z[787] = self.z[736] + self.z[276]*self.z[751] - self.z[273]*self.z[748] - self.z[279]*self.z[754] - self.z[284]*self.z[770]
        self.z[648] = id22*self.z[42]
        self.z[737] = self.z[47]*self.z[648]
        self.z[622] = ic22*self.z[42]
        self.z[618] = ic12*self.z[42]
        self.z[626] = ic23*self.z[42]
        self.z[663] = ie11*self.z[81]
        self.z[675] = ie12*self.z[82]
        self.z[688] = ie31*self.z[86]
        self.z[668] = ie12*self.z[81]
        self.z[679] = ie22*self.z[82]
        self.z[693] = ie23*self.z[86]
        self.z[713] = if22*self.z[82]
        self.z[673] = ie31*self.z[81]
        self.z[683] = ie23*self.z[82]
        self.z[698] = ie33*self.z[86]
        self.z[711] = if11*self.z[93]
        self.z[719] = if33*self.z[95]
        self.z[756] = self.z[40]*self.z[622] + self.z[40]*self.z[648] + self.z[43]*self.z[618] + self.z[45]*self.z[626] + self.z[78]*self.z[663] + self.z[78]*self.z[675] + self.z[78]*self.z[688] + self.z[84]*self.z[668] + self.z[84]*self.z[679] + self.z[84]*self.z[693] + self.z[84]*self.z[713] + self.z[85]*self.z[673] + self.z[85]*self.z[683] + self.z[85]*self.z[698] + self.z[90]*self.z[711] + self.z[94]*self.z[719] - massc*(self.z[5]*self.z[99]*self.z[104]-self.z[103]*self.z[104]-self.z[108]*self.z[109]-self.z[6]*self.z[99]*self.z[109]) - masse*(self.z[13]*self.z[99]*self.z[115]-self.z[113]*self.z[114]-self.z[113]*self.z[122]-self.z[114]*self.z[124]-self.z[115]*self.z[117]-self.z[118]*self.z[119]-self.z[122]*self.z[124]-self.z[6]*self.z[99]*self.z[114]-self.z[6]*self.z[99]*self.z[122]-self.z[9]*self.z[100]*self.z[119]-self.z[9]*self.z[110]*self.z[119]-self.z[10]*self.z[100]*self.z[115]-self.z[10]*self.z[110]*self.z[115]-self.z[15]*self.z[99]*self.z[119]) - massf*(self.z[13]*self.z[99]*self.z[125]-self.z[113]*self.z[114]-self.z[113]*self.z[132]-self.z[114]*self.z[134]-self.z[125]*self.z[127]-self.z[128]*self.z[129]-self.z[132]*self.z[134]-self.z[6]*self.z[99]*self.z[114]-self.z[6]*self.z[99]*self.z[132]-self.z[9]*self.z[100]*self.z[129]-self.z[9]*self.z[110]*self.z[129]-self.z[10]*self.z[100]*self.z[125]-self.z[10]*self.z[110]*self.z[125]-self.z[15]*self.z[99]*self.z[129])
        self.z[768] = self.z[89]*self.z[713]
        self.z[788] = self.z[737] + self.z[276]*self.z[752] - self.z[273]*self.z[749] - self.z[279]*self.z[756] - self.z[284]*self.z[768]
        self.z[792] = self.z[747]*self.z[276] - self.z[279]*self.z[751]
        self.z[793] = -self.z[279]*self.z[759] - self.z[284]*self.z[771]
        self.z[797] = self.z[791] + self.z[275]*self.z[790] + self.z[281]*self.z[787] + self.z[283]*self.z[788] - self.z[278]*self.z[792] - self.z[286]*self.z[793]
        self.z[730] = self.z[80]*self.z[661] + self.z[80]*self.z[690] + self.z[83]*self.z[666] + self.z[83]*self.z[695] + self.z[87]*self.z[671] + self.z[87]*self.z[700] + self.z[92]*self.z[709] + self.z[96]*self.z[721] + masse*self.z[121]*(self.z[120]+self.z[9]*self.z[101]+self.z[9]*self.z[111]+self.z[15]*self.z[98]) + massf*self.z[131]*(self.z[130]+self.z[9]*self.z[101]+self.z[9]*self.z[111]+self.z[15]*self.z[98])
        self.z[765] = self.z[81]*self.z[661] + self.z[81]*self.z[690] + self.z[82]*self.z[666] + self.z[82]*self.z[695] + self.z[86]*self.z[671] + self.z[86]*self.z[700] + self.z[93]*self.z[709] + self.z[95]*self.z[721] + masse*self.z[119]*self.z[121] +massf*self.z[129]*self.z[131]
        self.z[777] = self.z[730] + self.z[274]*self.z[744] + self.z[280]*self.z[758] + self.z[282]*self.z[765] - self.z[277]*self.z[745]
        self.z[728] = self.z[39]*self.z[622] + self.z[39]*self.z[648] + self.z[44]*self.z[618] + self.z[46]*self.z[626] + self.z[80]*self.z[663] + self.z[80]*self.z[675] + self.z[80]*self.z[688] + self.z[83]*self.z[668] + self.z[83]*self.z[679] + self.z[83]*self.z[693] + self.z[83]*self.z[713] + self.z[87]*self.z[673] + self.z[87]*self.z[683] + self.z[87]*self.z[698] + self.z[92]*self.z[711] + self.z[96]*self.z[719] - massc*(self.z[5]*self.z[98]*self.z[104]-self.z[102]*self.z[104]-self.z[107]*self.z[109]-self.z[6]*self.z[98]*self.z[109]) - masse*(self.z[13]*self.z[98]*self.z[115]-self.z[112]*self.z[114]-self.z[112]*self.z[122]-self.z[114]*self.z[123]-self.z[115]*self.z[116]-self.z[119]*self.z[120]-self.z[122]*self.z[123]-self.z[6]*self.z[98]*self.z[114]-self.z[6]*self.z[98]*self.z[122]-self.z[9]*self.z[101]*self.z[119]-self.z[9]*self.z[111]*self.z[119]-self.z[10]*self.z[101]*self.z[115]-self.z[10]*self.z[111]*self.z[115]-self.z[15]*self.z[98]*self.z[119]) - massf*(self.z[13]*self.z[98]*self.z[125]-self.z[112]*self.z[114]-self.z[112]*self.z[132]-self.z[114]*self.z[133]-self.z[125]*self.z[126]-self.z[129]*self.z[130]-self.z[132]*self.z[133]-self.z[6]*self.z[98]*self.z[114]-self.z[6]*self.z[98]*self.z[132]-self.z[9]*self.z[101]*self.z[129]-self.z[9]*self.z[111]*self.z[129]-self.z[10]*self.z[101]*self.z[125]-self.z[10]*self.z[111]*self.z[125]-self.z[15]*self.z[98]*self.z[129])
        self.z[763] = self.z[42]*self.z[622] + self.z[42]*self.z[648] + self.z[81]*self.z[663] + self.z[81]*self.z[675] + self.z[81]*self.z[688] + self.z[82]*self.z[668] + self.z[82]*self.z[679] + self.z[82]*self.z[693] + self.z[82]*self.z[713] + self.z[86]*self.z[673] + self.z[86]*self.z[683] + self.z[86]*self.z[698] + self.z[93]*self.z[711] + self.z[95]*self.z[719] + massc*(pow(self.z[104],2)+pow(self.z[109],2)) + masse*(pow(self.z[114],2)+pow(self.z[115],2)+pow(self.z[119],2)+pow(self.z[122],2)+2*self.z[114]*self.z[122]) + massf*(pow(self.z[114],2)+pow(self.z[125],2)+pow(self.z[129],2)+pow(self.z[132],2)+2*self.z[114]*self.z[132])
        self.z[775] = self.z[728] + self.z[274]*self.z[749] + self.z[280]*self.z[756] + self.z[282]*self.z[763] - self.z[277]*self.z[752] - self.z[285]*self.z[768]
        self.z[784] = self.z[777] + self.z[275]*self.z[780] + self.z[281]*self.z[773] + self.z[283]*self.z[775] - self.z[278]*self.z[779] - self.z[286]*self.z[778]
        self.z[738] = self.z[47]*self.z[649]
        self.z[789] = self.z[738] - self.z[279]*self.z[757]
        self.z[796] = self.z[789] + self.z[276]*self.z[792] - self.z[273]*self.z[790] - self.z[279]*self.z[787] - self.z[284]*self.z[793]
        self.z[817] = self.z[783]*self.z[797] - self.z[784]*self.z[796]
        self.z[292] = self.z[1]*self.z[3]*q4p - self.z[2]*self.z[4]*q3p
        self.z[328] = self.z[1]*self.z[6]*q3p + self.z[2]*self.z[5]*q5p + self.z[6]*self.z[20]*q5p - self.z[5]*self.z[292]
        self.z[315] = self.z[5]*self.z[9]*q7p - self.z[6]*self.z[10]*q5p
        self.z[311] = self.z[5]*self.z[10]*q5p + self.z[6]*self.z[9]*q7p
        self.z[291] = -self.z[1]*self.z[4]*q4p - self.z[2]*self.z[3]*q3p
        self.z[376] = self.z[1]*self.z[15]*q3p + self.z[10]*self.z[18]*q7p + self.z[2]*self.z[315] + self.z[16]*self.z[292] +self.z[20]*self.z[311] - self.z[9]*self.z[291]
        self.z[289] = self.z[1]*self.z[4]*q3p + self.z[2]*self.z[3]*q4p
        self.z[288] = self.z[1]*self.z[3]*q3p - self.z[2]*self.z[4]*q4p
        self.z[333] = self.z[2]*self.z[15]*q3p + self.z[10]*self.z[17]*q7p + self.z[16]*self.z[289] + self.z[19]*self.z[311] -self.z[1]*self.z[315] - self.z[9]*self.z[288]
        self.z[326] = self.z[1]*self.z[5]*q5p + self.z[5]*self.z[289] - self.z[2]*self.z[6]*q3p - self.z[6]*self.z[19]*q5p
        self.z[377] = self.z[24]*self.z[328] + self.z[27]*self.z[376] + self.z[28]*self.z[333] + self.z[62]*self.z[326]
        self.z[357] = self.z[9]*self.z[289] + self.z[16]*self.z[288] + self.z[17]*self.z[311] - self.z[10]*self.z[19]*q7p
        self.z[358] = self.z[10]*self.z[20]*q7p - self.z[9]*self.z[292] - self.z[16]*self.z[291] - self.z[18]*self.z[311]
        self.z[359] = -self.z[3]*self.z[10]*q7p - self.z[3]*self.z[16]*q4p - self.z[4]*self.z[9]*q4p - self.z[4]*self.z[311]
        self.z[330] = -self.z[3]*self.z[6]*q5p - self.z[4]*self.z[5]*q4p
        self.z[360] = self.z[27]*self.z[357] + self.z[28]*self.z[358] + self.z[29]*self.z[359] + self.z[64]*self.z[326] + self.z[66]*self.z[328] + self.z[68]*self.z[330]
        self.z[327] = self.z[1]*self.z[311] + self.z[15]*self.z[289] + self.z[19]*self.z[315] - self.z[2]*self.z[16]*q3p
        self.z[329] = self.z[1]*self.z[16]*q3p + self.z[2]*self.z[311] - self.z[15]*self.z[292] - self.z[20]*self.z[315]
        self.z[331] = self.z[3]*self.z[315] - self.z[4]*self.z[15]*q4p
        self.z[332] = self.z[27]*self.z[327] + self.z[28]*self.z[329] + self.z[29]*self.z[331] + self.z[65]*self.z[326] + self.z[67]*self.z[328] + self.z[70]*self.z[330]
        self.z[318] = self.z[5]*self.z[9]*q5p - self.z[6]*self.z[10]*q7p
        self.z[317] = -self.z[5]*self.z[10]*q7p - self.z[6]*self.z[9]*q5p
        self.z[389] = self.z[2]*self.z[13]*q3p + self.z[9]*self.z[17]*q7p + self.z[10]*self.z[288] + self.z[14]*self.z[289] +self.z[19]*self.z[318] - self.z[1]*self.z[317]
        self.z[378] = -self.z[1]*self.z[13]*q3p - self.z[9]*self.z[18]*q7p - self.z[2]*self.z[317] - self.z[10]*self.z[291] -self.z[14]*self.z[292] - self.z[20]*self.z[318]
        self.z[390] = self.z[3]*self.z[318] - self.z[3]*self.z[10]*q4p - self.z[4]*self.z[9]*q7p - self.z[4]*self.z[14]*q4p
        self.z[391] = self.z[27]*self.z[389] + self.z[28]*self.z[378] + self.z[29]*self.z[390] + self.z[55]*self.z[328] + self.z[63]*self.z[326] + self.z[69]*self.z[330]
        self.z[455] = q3p*self.z[377] + q4p*self.z[360] + q5p*self.z[332] + q7p*self.z[391]
        self.z[664] = ie11*self.z[455]
        self.z[374] = self.z[5]*self.z[292] - self.z[1]*self.z[6]*q3p - self.z[2]*self.z[5]*q5p - self.z[6]*self.z[20]*q5p
        self.z[321] = self.z[1]*self.z[13]*q3p + self.z[9]*self.z[18]*q7p + self.z[2]*self.z[317] + self.z[10]*self.z[291] +self.z[14]*self.z[292] + self.z[20]*self.z[318]
        self.z[319] = self.z[1]*self.z[317] - self.z[2]*self.z[13]*q3p - self.z[9]*self.z[17]*q7p - self.z[10]*self.z[288] -self.z[14]*self.z[289] - self.z[19]*self.z[318]
        self.z[375] = self.z[21]*self.z[374] + self.z[22]*self.z[326] + self.z[27]*self.z[321] + self.z[73]*self.z[319]
        self.z[353] = self.z[5]*self.z[288] - self.z[6]*self.z[17]*q5p
        self.z[354] = self.z[5]*self.z[291] - self.z[6]*self.z[18]*q5p
        self.z[355] = self.z[3]*self.z[5]*q4p - self.z[4]*self.z[6]*q5p
        self.z[323] = self.z[3]*self.z[10]*q4p + self.z[4]*self.z[9]*q7p + self.z[4]*self.z[14]*q4p - self.z[3]*self.z[318]
        self.z[356] = self.z[21]*self.z[353] + self.z[72]*self.z[319] - self.z[22]*self.z[354] - self.z[23]*self.z[355] - self.z[75]*self.z[321] - self.z[76]*self.z[323]
        self.z[320] = -self.z[1]*self.z[6]*q5p - self.z[2]*self.z[5]*q3p - self.z[5]*self.z[19]*q5p - self.z[6]*self.z[289]
        self.z[322] = self.z[1]*self.z[5]*q3p + self.z[5]*self.z[20]*q5p + self.z[6]*self.z[292] - self.z[2]*self.z[6]*q5p
        self.z[324] = self.z[3]*self.z[5]*q5p - self.z[4]*self.z[6]*q4p
        self.z[325] = self.z[21]*self.z[320] + self.z[22]*self.z[322] + self.z[71]*self.z[319] + self.z[74]*self.z[321] - self.z[23]*self.z[324] - self.z[77]*self.z[323]
        self.z[456] = q3p*self.z[375] + q4p*self.z[356] + q5p*self.z[325]
        self.z[678] = ie12*self.z[456]
        self.z[316] = self.z[9]*self.z[291] - self.z[1]*self.z[15]*q3p - self.z[10]*self.z[18]*q7p - self.z[2]*self.z[315] -self.z[16]*self.z[292] - self.z[20]*self.z[311]
        self.z[379] = self.z[21]*self.z[316] + self.z[24]*self.z[378] + self.z[25]*self.z[319] + self.z[55]*self.z[333]
        self.z[361] = self.z[9]*self.z[19]*q7p + self.z[10]*self.z[289] - self.z[14]*self.z[288] - self.z[17]*self.z[318]
        self.z[362] = self.z[14]*self.z[291] + self.z[18]*self.z[318] - self.z[9]*self.z[20]*q7p - self.z[10]*self.z[292]
        self.z[363] = self.z[3]*self.z[9]*q7p + self.z[3]*self.z[14]*q4p + self.z[4]*self.z[318] - self.z[4]*self.z[10]*q4p
        self.z[312] = self.z[3]*self.z[9]*q4p + self.z[3]*self.z[311] - self.z[4]*self.z[10]*q7p - self.z[4]*self.z[16]*q4p
        self.z[364] = self.z[24]*self.z[361] + self.z[25]*self.z[362] + self.z[26]*self.z[363] + self.z[57]*self.z[333] + self.z[58]*self.z[316] + self.z[60]*self.z[312]
        self.z[334] = self.z[2]*self.z[14]*q3p - self.z[1]*self.z[318] - self.z[13]*self.z[289] - self.z[19]*self.z[317]
        self.z[335] = self.z[13]*self.z[292] + self.z[20]*self.z[317] - self.z[1]*self.z[14]*q3p - self.z[2]*self.z[318]
        self.z[336] = self.z[3]*self.z[317] - self.z[4]*self.z[13]*q4p
        self.z[337] = self.z[24]*self.z[334] + self.z[25]*self.z[335] + self.z[56]*self.z[333] + self.z[59]*self.z[316] - self.z[26]*self.z[336] - self.z[61]*self.z[312]
        self.z[392] = 2*self.z[24]*self.z[333] + 2*self.z[25]*self.z[316] + 2*self.z[26]*self.z[312]
        self.z[457] = q3p*self.z[379] + q4p*self.z[364] + q5p*self.z[337] + q7p*self.z[392]
        self.z[691] = ie31*self.z[457]
        self.z[453] = self.z[82]*q5p + self.z[83]*q4p + self.z[84]*q3p
        self.z[452] = self.z[78]*q3p + self.z[79]*q7p + self.z[80]*q4p + self.z[81]*q5p
        self.z[454] = self.z[85]*q3p + self.z[86]*q5p + self.z[87]*q4p + self.z[88]*q7p
        self.z[659] = ie23*self.z[453] + ie31*self.z[452] + ie33*self.z[454]
        self.z[658] = ie12*self.z[452] + ie22*self.z[453] + ie23*self.z[454]
        self.z[704] = self.z[453]*self.z[659] - self.z[454]*self.z[658]
        self.z[674] = ie31*self.z[455]
        self.z[686] = ie23*self.z[456]
        self.z[701] = ie33*self.z[457]
        self.z[657] = ie11*self.z[452] + ie12*self.z[453] + ie31*self.z[454]
        self.z[702] = self.z[452]*self.z[658] - self.z[453]*self.z[657]
        self.z[380] = self.z[11]*self.z[377] - self.z[11]*self.z[85]*q8p - self.z[12]*self.z[78]*q8p - self.z[12]*self.z[379]
        self.z[365] = self.z[11]*self.z[360] - self.z[11]*self.z[87]*q8p - self.z[12]*self.z[80]*q8p - self.z[12]*self.z[364]
        self.z[338] = self.z[11]*self.z[332] - self.z[11]*self.z[86]*q8p - self.z[12]*self.z[81]*q8p - self.z[12]*self.z[337]
        self.z[393] = self.z[11]*self.z[391] - self.z[11]*self.z[88]*q8p - self.z[12]*self.z[79]*q8p - self.z[12]*self.z[392]
        self.z[461] = q3p*self.z[380] + q4p*self.z[365] + q5p*self.z[338] + q7p*self.z[393]
        self.z[712] = if11*self.z[461]
        self.z[459] = self.z[82]*q5p + self.z[83]*q4p + self.z[84]*q3p + self.z[89]*q8p
        self.z[460] = self.z[94]*q3p + self.z[95]*q5p + self.z[96]*q4p + self.z[97]*q7p
        self.z[707] = if33*self.z[460]
        self.z[706] = if22*self.z[459]
        self.z[725] = self.z[459]*self.z[707] - self.z[460]*self.z[706]
        self.z[381] = self.z[11]*self.z[78]*q8p + self.z[11]*self.z[379] + self.z[12]*self.z[377] - self.z[12]*self.z[85]*q8p
        self.z[366] = self.z[11]*self.z[80]*q8p + self.z[11]*self.z[364] + self.z[12]*self.z[360] - self.z[12]*self.z[87]*q8p
        self.z[340] = self.z[11]*self.z[81]*q8p + self.z[11]*self.z[337] + self.z[12]*self.z[332] - self.z[12]*self.z[86]*q8p
        self.z[394] = self.z[11]*self.z[79]*q8p + self.z[11]*self.z[392] + self.z[12]*self.z[391] - self.z[12]*self.z[88]*q8p
        self.z[462] = q3p*self.z[381] + q4p*self.z[366] + q5p*self.z[340] + q7p*self.z[394]
        self.z[722] = if33*self.z[462]
        self.z[458] = self.z[90]*q3p + self.z[91]*q7p + self.z[92]*q4p + self.z[93]*q5p
        self.z[705] = if11*self.z[458]
        self.z[723] = self.z[458]*self.z[706] - self.z[459]*self.z[705]
        self.z[440] = self.z[39]*q4p + self.z[40]*q3p + self.z[42]*q5p
        self.z[485] = -self.z[112]*q4p - self.z[113]*q3p - self.z[114]*q5p
        self.z[441] = self.z[45]*q3p + self.z[46]*q4p
        self.z[484] = self.z[110]*q3p + self.z[111]*q4p
        self.z[498] = self.z[440]*self.z[485] - self.z[441]*self.z[484]
        self.z[384] = self.z[1]*self.z[19]*q3p + self.z[1]*self.z[292] + self.z[2]*self.z[289] - self.z[2]*self.z[20]*q3p
        self.z[369] = self.z[1]*self.z[288] - self.z[1]*self.z[18]*q3p - self.z[2]*self.z[17]*q3p - self.z[2]*self.z[291]
        self.z[468] = rR*(q3p*self.z[384]+q4p*self.z[369])
        self.z[464] = self.z[100]*q3p + self.z[101]*q4p
        self.z[471] = self.z[41]*q3p
        self.z[478] = -self.z[468] - self.z[464]*self.z[471]
        self.z[492] = l3*self.z[379] - l4*self.z[377]
        self.z[494] = l3*self.z[364] - l4*self.z[360]
        self.z[493] = l3*self.z[337] - l4*self.z[332]
        self.z[495] = l3*self.z[392] - l4*self.z[391]
        self.z[496] = q3p*self.z[492] + q4p*self.z[494] + q5p*self.z[493] + q7p*self.z[495]
        self.z[486] = self.z[115]*q5p + self.z[116]*q4p + self.z[117]*q3p
        self.z[488] = -self.z[122]*q5p - self.z[123]*q4p - self.z[124]*q3p
        self.z[502] = self.z[496] + self.z[454]*self.z[486] - self.z[452]*self.z[488]
        self.z[294] = self.z[17]*self.z[292] + self.z[20]*self.z[288] - self.z[18]*self.z[289] - self.z[19]*self.z[291]
        self.z[304] = 2*self.z[19]*self.z[289] + 2*self.z[20]*self.z[292] - 2*self.z[3]*self.z[4]*q4p
        self.z[472] = rR*(q3p*self.z[294]+q4p*self.z[304])
        self.z[463] = -self.z[98]*q4p - self.z[99]*q3p
        self.z[479] = self.z[472] + self.z[463]*self.z[471]
        self.z[296] = self.z[1]*self.z[17]*q3p + self.z[1]*self.z[291] + self.z[2]*self.z[288] - self.z[2]*self.z[18]*q3p
        self.z[405] = self.z[5]*self.z[37]*q5p + self.z[5]*self.z[296] + self.z[6]*self.z[294] - self.z[6]*self.z[41]*q5p
        self.z[412] = self.z[5]*self.z[38]*q5p + self.z[6]*self.z[304]
        self.z[489] = d1*(q3p*self.z[405]+q4p*self.z[412])
        self.z[439] = self.z[43]*q3p + self.z[44]*q4p
        self.z[499] = self.z[489] - self.z[439]*self.z[485]
        self.z[469] = self.z[37]*q3p + self.z[38]*q4p
        self.z[470] = self.z[39]*q4p + self.z[40]*q3p
        self.z[480] = self.z[464]*self.z[469] - self.z[463]*self.z[470]
        self.z[386] = d3*self.z[379] - d2*self.z[377]
        self.z[371] = d3*self.z[364] - d2*self.z[360]
        self.z[349] = d3*self.z[337] - d2*self.z[332]
        self.z[398] = d3*self.z[392] - d2*self.z[391]
        self.z[508] = q3p*self.z[386] + q4p*self.z[371] + q5p*self.z[349] + q7p*self.z[398]
        self.z[504] = self.z[125]*q5p + self.z[126]*q4p + self.z[127]*q3p
        self.z[506] = -self.z[132]*q5p - self.z[133]*q4p - self.z[134]*q3p
        self.z[511] = self.z[508] + self.z[454]*self.z[504] - self.z[452]*self.z[506]
        self.z[746] = self.z[79]*self.z[664] + self.z[79]*self.z[678] + self.z[79]*self.z[691] + self.z[79]*self.z[704] + self.z[88]*self.z[674] + self.z[88]*self.z[686] + self.z[88]*self.z[701] + self.z[88]*self.z[702] + self.z[91]*self.z[712] + self.z[91]*self.z[725] + self.z[97]*self.z[722] + self.z[97]*self.z[723] - masse*self.z[121]*(self.z[10]*self.z[498]+self.z[15]*self.z[478]-self.z[502]-self.z[9]*self.z[479]-self.z[9]*self.z[499]-self.z[16]*self.z[480]) - massf*self.z[131]*(self.z[10]*self.z[498]+self.z[15]*self.z[478]-self.z[511]-self.z[9]*self.z[479]-self.z[9]*self.z[499]-self.z[16]*self.z[480])
        self.z[442] = self.z[5]*self.z[294] - self.z[5]*self.z[41]*q5p - self.z[6]*self.z[37]*q5p - self.z[6]*self.z[296]
        self.z[474] = l1*self.z[405] - l2*self.z[442]
        self.z[443] = self.z[5]*self.z[304] - self.z[6]*self.z[38]*q5p
        self.z[475] = l1*self.z[412] - l2*self.z[443]
        self.z[476] = q3p*self.z[474] + q4p*self.z[475]
        self.z[465] = self.z[102]*q4p + self.z[103]*q3p + self.z[104]*q5p
        self.z[467] = -self.z[107]*q4p - self.z[108]*q3p - self.z[109]*q5p
        self.z[482] = self.z[476] + self.z[441]*self.z[465] - self.z[439]*self.z[467]
        self.z[466] = self.z[105]*q3p + self.z[106]*q4p
        self.z[477] = l1*(q3p*self.z[384]+q4p*self.z[369])
        self.z[483] = self.z[439]*self.z[466] - self.z[477] - self.z[440]*self.z[465]
        self.z[473] = l2*(q3p*self.z[384]+q4p*self.z[369])
        self.z[481] = self.z[473] + self.z[440]*self.z[467] - self.z[441]*self.z[466]
        self.z[491] = l4*(q3p*self.z[375]+q4p*self.z[356]+q5p*self.z[325])
        self.z[487] = self.z[118]*q3p + self.z[119]*q5p + self.z[120]*q4p + self.z[121]*q7p
        self.z[501] = self.z[491] + self.z[453]*self.z[488] - self.z[454]*self.z[487]
        self.z[490] = d1*(q3p*self.z[384]+q4p*self.z[369])
        self.z[500] = self.z[439]*self.z[484] - self.z[490]
        self.z[497] = l3*(q3p*self.z[375]+q4p*self.z[356]+q5p*self.z[325])
        self.z[503] = self.z[452]*self.z[487] - self.z[497] - self.z[453]*self.z[486]
        self.z[507] = d2*(q3p*self.z[375]+q4p*self.z[356]+q5p*self.z[325])
        self.z[505] = self.z[128]*q3p + self.z[129]*q5p + self.z[130]*q4p + self.z[131]*q7p
        self.z[510] = self.z[507] + self.z[453]*self.z[506] - self.z[454]*self.z[505]
        self.z[509] = d3*(q3p*self.z[375]+q4p*self.z[356]+q5p*self.z[325])
        self.z[512] = self.z[452]*self.z[505] - self.z[509] - self.z[453]*self.z[504]
        self.z[750] = -massd*(self.z[17]*self.z[479]-self.z[1]*self.z[478]-self.z[19]*self.z[480]) - massc*(self.z[17]*self.z[479]+self.z[17]*self.z[482]-self.z[1]*self.z[478]-self.z[19]*self.z[480]-self.z[27]*self.z[483]-self.z[71]*self.z[481]) -masse*(self.z[17]*self.z[479]+self.z[17]*self.z[499]-self.z[1]*self.z[478]-self.z[19]*self.z[480]-self.z[21]*self.z[501]-self.z[24]*self.z[502]-self.z[27]*self.z[500]-self.z[27]*self.z[503]-self.z[71]*self.z[498]) - massf*(self.z[17]*self.z[479]+self.z[17]*self.z[499]-self.z[1]*self.z[478]-self.z[19]*self.z[480]-self.z[21]*self.z[510]-self.z[24]*self.z[511]-self.z[27]*self.z[500]-self.z[27]*self.z[512]-self.z[71]*self.z[498])
        self.z[444] = q3p*self.z[442] + q4p*self.z[443]
        self.z[612] = ic12*self.z[444]
        self.z[445] = q3p*self.z[384] + q4p*self.z[369]
        self.z[623] = ic22*self.z[445]
        self.z[446] = q3p*self.z[405] + q4p*self.z[412]
        self.z[633] = ic23*self.z[446]
        self.z[604] = ic11*self.z[439] + ic12*self.z[440] + ic31*self.z[441]
        self.z[606] = ic23*self.z[440] + ic31*self.z[439] + ic33*self.z[441]
        self.z[638] = self.z[441]*self.z[604] - self.z[439]*self.z[606]
        self.z[650] = id22*self.z[445]
        self.z[449] = self.z[53]*q3p + self.z[54]*q4p
        self.z[447] = self.z[51]*q3p + self.z[52]*q4p
        self.z[640] = id11*self.z[447]
        self.z[642] = id33*self.z[449]
        self.z[655] = self.z[449]*self.z[640] - self.z[447]*self.z[642]
        self.z[609] = ic11*self.z[444]
        self.z[619] = ic12*self.z[445]
        self.z[630] = ic31*self.z[446]
        self.z[605] = ic12*self.z[439] + ic22*self.z[440] + ic23*self.z[441]
        self.z[639] = self.z[440]*self.z[606] - self.z[441]*self.z[605]
        self.z[615] = ic31*self.z[444]
        self.z[627] = ic23*self.z[445]
        self.z[636] = ic33*self.z[446]
        self.z[637] = self.z[439]*self.z[605] - self.z[440]*self.z[604]
        self.z[295] = -self.z[5]*self.z[8]*q5p - self.z[5]*self.z[8]*q6p - self.z[6]*self.z[7]*q5p - self.z[6]*self.z[7]*q6p
        self.z[297] = self.z[6]*self.z[8]*q5p + self.z[6]*self.z[8]*q6p - self.z[5]*self.z[7]*q5p - self.z[5]*self.z[7]*q6p
        self.z[298] = self.z[37]*self.z[295] + self.z[41]*self.z[297] + self.z[48]*self.z[294] + self.z[50]*self.z[296]
        self.z[305] = self.z[38]*self.z[295] + self.z[48]*self.z[304]
        self.z[450] = q3p*self.z[298] + q4p*self.z[305]
        self.z[645] = id11*self.z[450]
        self.z[448] = self.z[39]*q4p + self.z[40]*q3p + self.z[42]*q5p + self.z[47]*q6p
        self.z[641] = id22*self.z[448]
        self.z[656] = self.z[448]*self.z[642] - self.z[449]*self.z[641]
        self.z[299] = self.z[5]*self.z[7]*q5p + self.z[5]*self.z[7]*q6p - self.z[6]*self.z[8]*q5p - self.z[6]*self.z[8]*q6p
        self.z[300] = self.z[37]*self.z[299] + self.z[41]*self.z[295] + self.z[48]*self.z[296] + self.z[49]*self.z[294]
        self.z[306] = self.z[38]*self.z[299] + self.z[49]*self.z[304]
        self.z[451] = q3p*self.z[300] + q4p*self.z[306]
        self.z[653] = id33*self.z[451]
        self.z[654] = self.z[447]*self.z[641] - self.z[448]*self.z[640]
        self.z[669] = ie12*self.z[455]
        self.z[682] = ie22*self.z[456]
        self.z[696] = ie23*self.z[457]
        self.z[703] = self.z[454]*self.z[657] - self.z[452]*self.z[659]
        self.z[717] = if22*self.z[456]
        self.z[724] = self.z[460]*self.z[705] - self.z[458]*self.z[707]
        self.z[760] = self.z[40]*self.z[612] + self.z[40]*self.z[623] + self.z[40]*self.z[633] + self.z[40]*self.z[638] + self.z[40]*self.z[650] + self.z[40]*self.z[655] + self.z[43]*self.z[609] + self.z[43]*self.z[619] + self.z[43]*self.z[630] + self.z[43]*self.z[639] + self.z[45]*self.z[615] + self.z[45]*self.z[627] + self.z[45]*self.z[636] + self.z[45]*self.z[637] + self.z[51]*self.z[645] + self.z[51]*self.z[656] + self.z[53]*self.z[653] + self.z[53]*self.z[654] + self.z[78]*self.z[664] + self.z[78]*self.z[678] + self.z[78]*self.z[691] + self.z[78]*self.z[704] + self.z[84]*self.z[669] + self.z[84]*self.z[682] + self.z[84]*self.z[696] + self.z[84]*self.z[703] + self.z[84]*self.z[717] + self.z[84]*self.z[724] + self.z[85]*self.z[674] + self.z[85]*self.z[686] + self.z[85]*self.z[701] + self.z[85]*self.z[702] + self.z[90]*self.z[712] + self.z[90]*self.z[725] + self.z[94]*self.z[722] + self.z[94]*self.z[723] + massc*(self.z[100]*self.z[479]+self.z[100]*self.z[482]+self.z[103]*self.z[481]+self.z[105]*self.z[479]+self.z[105]*self.z[482]+self.z[5]*self.z[103]*self.z[478]-self.z[99]*self.z[478]-self.z[108]*self.z[483]-self.z[5]*self.z[99]*self.z[481]-self.z[5]*self.z[108]*self.z[480]-self.z[6]*self.z[99]*self.z[483]-self.z[6]*self.z[103]*self.z[480]-self.z[6]*self.z[108]*self.z[478]) - massd*(self.z[99]*self.z[478]-self.z[100]*self.z[479]) - masse*(self.z[99]*self.z[478]+self.z[113]*self.z[500]+self.z[113]*self.z[503]+self.z[124]*self.z[500]+self.z[124]*self.z[503]+self.z[5]*self.z[99]*self.z[498]+self.z[5]*self.z[113]*self.z[480]+self.z[5]*self.z[124]*self.z[480]+self.z[6]*self.z[99]*self.z[500]+self.z[6]*self.z[99]*self.z[503]+self.z[6]*self.z[113]*self.z[478]+self.z[6]*self.z[124]*self.z[478]+self.z[10]*self.z[118]*self.z[498]+self.z[13]*self.z[99]*self.z[501]+self.z[14]*self.z[117]*self.z[480]+self.z[15]*self.z[118]*self.z[478]-self.z[100]*self.z[479]-self.z[100]*self.z[499]-self.z[110]*self.z[479]-self.z[110]*self.z[499]-self.z[117]*self.z[501]-self.z[118]*self.z[502]-self.z[9]*self.z[100]*self.z[502]-self.z[9]*self.z[110]*self.z[502]-self.z[9]*self.z[117]*self.z[498]-self.z[9]*self.z[118]*self.z[479]-self.z[9]*self.z[118]*self.z[499]-self.z[10]*self.z[100]*self.z[501]-self.z[10]*self.z[110]*self.z[501]-self.z[10]*self.z[117]*self.z[479]-self.z[10]*self.z[117]*self.z[499]-self.z[13]*self.z[117]*self.z[478]-self.z[15]*self.z[99]*self.z[502]-self.z[16]*self.z[118]*self.z[480]) -massf*(self.z[99]*self.z[478]+self.z[113]*self.z[500]+self.z[113]*self.z[512]+self.z[134]*self.z[500]+self.z[134]*self.z[512]+self.z[5]*self.z[99]*self.z[498]+self.z[5]*self.z[113]*self.z[480]+self.z[5]*self.z[134]*self.z[480]+self.z[6]*self.z[99]*self.z[500]+self.z[6]*self.z[99]*self.z[512]+self.z[6]*self.z[113]*self.z[478]+self.z[6]*self.z[134]*self.z[478]+self.z[10]*self.z[128]*self.z[498]+self.z[13]*self.z[99]*self.z[510]+self.z[14]*self.z[127]*self.z[480]+self.z[15]*self.z[128]*self.z[478]-self.z[100]*self.z[479]-self.z[100]*self.z[499]-self.z[110]*self.z[479]-self.z[110]*self.z[499]-self.z[127]*self.z[510]-self.z[128]*self.z[511]-self.z[9]*self.z[100]*self.z[511]-self.z[9]*self.z[110]*self.z[511]-self.z[9]*self.z[127]*self.z[498]-self.z[9]*self.z[128]*self.z[479]-self.z[9]*self.z[128]*self.z[499]-self.z[10]*self.z[100]*self.z[510]-self.z[10]*self.z[110]*self.z[510]-self.z[10]*self.z[127]*self.z[479]-self.z[10]*self.z[127]*self.z[499]-self.z[13]*self.z[127]*self.z[478]-self.z[15]*self.z[99]*self.z[511]-self.z[16]*self.z[128]*self.z[480])
        self.z[767] = self.z[81]*self.z[664] + self.z[81]*self.z[678] + self.z[81]*self.z[691] + self.z[81]*self.z[704] + self.z[82]*self.z[669] + self.z[82]*self.z[682] + self.z[82]*self.z[696] + self.z[82]*self.z[703] + self.z[82]*self.z[717] + self.z[82]*self.z[724] + self.z[86]*self.z[674] + self.z[86]*self.z[686] + self.z[86]*self.z[701] + self.z[86]*self.z[702] + self.z[93]*self.z[712] + self.z[93]*self.z[725] + self.z[95]*self.z[722] + self.z[95]*self.z[723] + self.z[42]*(self.z[650]+self.z[655]) +self.z[42]*(self.z[612]+self.z[623]+self.z[633]+self.z[638]) + massc*(self.z[104]*self.z[481]+self.z[5]*self.z[104]*self.z[478]-self.z[109]*self.z[483]-self.z[5]*self.z[109]*self.z[480]-self.z[6]*self.z[104]*self.z[480]-self.z[6]*self.z[109]*self.z[478]) - masse*(self.z[114]*self.z[500]+self.z[114]*self.z[503]+self.z[122]*self.z[500]+self.z[122]*self.z[503]+self.z[5]*self.z[114]*self.z[480]+self.z[5]*self.z[122]*self.z[480]+self.z[6]*self.z[114]*self.z[478]+self.z[6]*self.z[122]*self.z[478]+self.z[10]*self.z[119]*self.z[498]+self.z[14]*self.z[115]*self.z[480]+self.z[15]*self.z[119]*self.z[478]-self.z[115]*self.z[501]-self.z[119]*self.z[502]-self.z[9]*self.z[115]*self.z[498]-self.z[9]*self.z[119]*self.z[479]-self.z[9]*self.z[119]*self.z[499]-self.z[10]*self.z[115]*self.z[479]-self.z[10]*self.z[115]*self.z[499]-self.z[13]*self.z[115]*self.z[478]-self.z[16]*self.z[119]*self.z[480]) - massf*(self.z[114]*self.z[500]+self.z[114]*self.z[512]+self.z[132]*self.z[500]+self.z[132]*self.z[512]+self.z[5]*self.z[114]*self.z[480]+self.z[5]*self.z[132]*self.z[480]+self.z[6]*self.z[114]*self.z[478]+self.z[6]*self.z[132]*self.z[478]+self.z[10]*self.z[129]*self.z[498]+self.z[14]*self.z[125]*self.z[480]+self.z[15]*self.z[129]*self.z[478]-self.z[125]*self.z[510]-self.z[129]*self.z[511]-self.z[9]*self.z[125]*self.z[498]-self.z[9]*self.z[129]*self.z[479]-self.z[9]*self.z[129]*self.z[499]-self.z[10]*self.z[125]*self.z[479]-self.z[10]*self.z[125]*self.z[499]-self.z[13]*self.z[125]*self.z[478]-self.z[16]*self.z[129]*self.z[480])
        self.z[753] = massd*(self.z[2]*self.z[478]+self.z[18]*self.z[479]-self.z[20]*self.z[480]) - massc*(self.z[20]*self.z[480]-self.z[2]*self.z[478]-self.z[18]*self.z[479]-self.z[18]*self.z[482]-self.z[28]*self.z[483]-self.z[74]*self.z[481]) -masse*(self.z[20]*self.z[480]-self.z[2]*self.z[478]-self.z[18]*self.z[479]-self.z[18]*self.z[499]-self.z[22]*self.z[501]-self.z[25]*self.z[502]-self.z[28]*self.z[500]-self.z[28]*self.z[503]-self.z[74]*self.z[498]) - massf*(self.z[20]*self.z[480]-self.z[2]*self.z[478]-self.z[18]*self.z[479]-self.z[18]*self.z[499]-self.z[22]*self.z[510]-self.z[25]*self.z[511]-self.z[28]*self.z[500]-self.z[28]*self.z[512]-self.z[74]*self.z[498])
        self.z[772] = self.z[89]*(self.z[717]+self.z[724])
        self.z[807] = self.z[746] + self.z[275]*self.z[750] + self.z[281]*self.z[760] + self.z[283]*self.z[767] - self.z[278]*self.z[753] - self.z[286]*self.z[772]
        self.z[740] = self.z[79]*self.z[663] + self.z[79]*self.z[675] + self.z[79]*self.z[688] + self.z[88]*self.z[673] + self.z[88]*self.z[683] + self.z[88]*self.z[698] + self.z[91]*self.z[711] + self.z[97]*self.z[719] + masse*self.z[119]*self.z[121] +massf*self.z[129]*self.z[131]
        self.z[799] = self.z[740] + self.z[275]*self.z[749] + self.z[281]*self.z[756] + self.z[283]*self.z[763] - self.z[278]*self.z[752] - self.z[286]*self.z[768]
        self.z[418] = self.z[26]*self.z[359] + self.z[68]*self.z[312]
        self.z[419] = (self.z[30]*self.z[418]+self.z[26]*self.z[184]*self.z[312])/pow(self.z[30],1.5)
        self.z[313] = self.z[26]*self.z[312]/pow(self.z[30],0.5)
        self.z[420] = (self.z[31]*self.z[419]+2*self.z[187]*self.z[313])/pow(self.z[31],3)
        self.z[406] = self.z[3]*self.z[6]*q4p + self.z[4]*self.z[5]*q5p
        self.z[421] = ((self.z[31]*self.z[68]+2*self.z[26]*self.z[187])*self.z[313]+self.z[31]*(self.z[26]*self.z[419]+self.z[31]*self.z[359]+self.z[187]*self.z[312]))/pow(self.z[31],3)
        self.z[343] = (self.z[26]*self.z[313]+self.z[31]*self.z[312])/pow(self.z[31],2)
        self.z[422] = rR*self.z[3]*q4p + rF*self.z[420] + d1*self.z[406] + d3*self.z[363] - d2*self.z[355] -self.z[34]*self.z[359] - self.z[199]*self.z[312] - rF*self.z[26]*self.z[421] - rF*self.z[68]*self.z[343]
        self.z[428] = self.z[26]*self.z[331] + self.z[70]*self.z[312]
        self.z[429] = (self.z[30]*self.z[428]+self.z[26]*self.z[186]*self.z[312])/pow(self.z[30],1.5)
        self.z[430] = (self.z[31]*self.z[429]+2*self.z[189]*self.z[313])/pow(self.z[31],3)
        self.z[431] = ((self.z[31]*self.z[70]+2*self.z[26]*self.z[189])*self.z[313]+self.z[31]*(self.z[26]*self.z[429]+self.z[31]*self.z[331]+self.z[189]*self.z[312]))/pow(self.z[31],3)
        self.z[432] = rF*self.z[430] - d1*self.z[330] - d2*self.z[324] - d3*self.z[336] - self.z[34]*self.z[331] -self.z[201]*self.z[312] - rF*self.z[26]*self.z[431] - rF*self.z[70]*self.z[343]
        self.z[423] = self.z[26]*self.z[390] + self.z[69]*self.z[312]
        self.z[424] = (self.z[30]*self.z[423]+self.z[26]*self.z[185]*self.z[312])/pow(self.z[30],1.5)
        self.z[425] = (self.z[31]*self.z[424]+2*self.z[188]*self.z[313])/pow(self.z[31],3)
        self.z[426] = ((self.z[31]*self.z[69]+2*self.z[26]*self.z[188])*self.z[313]+self.z[31]*(self.z[26]*self.z[424]+self.z[31]*self.z[390]+self.z[188]*self.z[312]))/pow(self.z[31],3)
        self.z[427] = rF*self.z[425] + d3*self.z[312] - self.z[34]*self.z[390] - self.z[200]*self.z[312] - rF*self.z[26]*self.z[426] - rF*self.z[69]*self.z[343]
        self.z[433] = u4*self.z[422] + u5*self.z[432] + u7*self.z[427]
        self.z[437] = self.z[271]*self.z[433]/self.z[234]
        self.z[742] = self.z[79]*self.z[660] + self.z[79]*self.z[677] + self.z[79]*self.z[687] + self.z[88]*self.z[670] + self.z[88]*self.z[685] + self.z[88]*self.z[697] + self.z[91]*self.z[708] + self.z[97]*self.z[718] + masse*self.z[121]*(self.z[118]+self.z[9]*self.z[100]+self.z[9]*self.z[110]+self.z[15]*self.z[99]) + massf*self.z[131]*(self.z[128]+self.z[9]*self.z[100]+self.z[9]*self.z[110]+self.z[15]*self.z[99])
        self.z[801] = self.z[742] + self.z[275]*self.z[748] + self.z[281]*self.z[754] + self.z[283]*self.z[761] - self.z[278]*self.z[751] - self.z[286]*self.z[770]
        self.z[290] = self.z[4]*self.z[17]*q4p - self.z[3]*self.z[19]*q4p - self.z[3]*self.z[288] - self.z[4]*self.z[289]
        self.z[293] = self.z[3]*self.z[20]*q4p + self.z[3]*self.z[291] + self.z[4]*self.z[292] - self.z[4]*self.z[18]*q4p
        self.z[301] = rR*(self.z[48]*self.z[298]+self.z[49]*self.z[300]+self.z[51]*self.z[295]+self.z[53]*self.z[299])
        self.z[302] = rR*self.z[294] - self.z[301]
        self.z[303] = self.z[3]*self.z[302] - self.z[4]*self.z[139]*q4p
        self.z[307] = rR*(self.z[48]*self.z[305]+self.z[49]*self.z[306]+self.z[52]*self.z[295]+self.z[54]*self.z[299])
        self.z[308] = rR*self.z[304] - self.z[307]
        self.z[309] = self.z[3]*self.z[308] - self.z[4]*self.z[140]*q4p
        self.z[310] = u1*self.z[290] + u2*self.z[293] + u3*self.z[303] + u4*self.z[309]
        self.z[314] = rF*self.z[313]/pow(self.z[31],2)
        self.z[339] = self.z[11]*self.z[321] - self.z[11]*self.z[28]*q8p - self.z[12]*self.z[22]*q8p - self.z[12]*self.z[328]
        self.z[341] = self.z[11]*self.z[22]*q8p + self.z[11]*self.z[328] + self.z[12]*self.z[321] - self.z[12]*self.z[28]*q8p
        self.z[382] = (self.z[25]*self.z[84]+self.z[90]*self.z[142]+self.z[94]*self.z[145])*self.z[314] + self.z[35]*(self.z[25]*self.z[375]+self.z[84]*self.z[316]+self.z[90]*self.z[339]+self.z[94]*self.z[341]+self.z[142]*self.z[380]+self.z[145]*self.z[381])
        self.z[345] = self.z[11]*self.z[317] - self.z[5]*self.z[12]*q5p - self.z[6]*self.z[11]*q8p - self.z[12]*self.z[13]*q8p
        self.z[383] = self.z[34]*self.z[381] + rF*self.z[94]*self.z[343]
        self.z[346] = self.z[11]*self.z[319] - self.z[11]*self.z[27]*q8p - self.z[12]*self.z[21]*q8p - self.z[12]*self.z[326]
        self.z[347] = self.z[11]*self.z[21]*q8p + self.z[11]*self.z[326] + self.z[12]*self.z[319] - self.z[12]*self.z[27]*q8p
        self.z[385] = (self.z[24]*self.z[84]+self.z[90]*self.z[141]+self.z[94]*self.z[144])*self.z[314] + self.z[35]*(self.z[24]*self.z[375]+self.z[84]*self.z[333]+self.z[90]*self.z[346]+self.z[94]*self.z[347]+self.z[141]*self.z[380]+self.z[144]*self.z[381])
        self.z[351] = self.z[5]*self.z[11]*q5p + self.z[11]*self.z[13]*q8p + self.z[12]*self.z[317] - self.z[6]*self.z[12]*q8p
        self.z[387] = self.z[34]*self.z[380] + rF*self.z[90]*self.z[343]
        self.z[388] = self.z[1]*self.z[382] + self.z[127]*self.z[317] + self.z[161]*self.z[345] + self.z[178]*self.z[383] + d2*self.z[13]*self.z[375] - self.z[1]*self.z[153]*q3p - self.z[2]*self.z[158]*q3p - self.z[5]*self.z[113]*q5p - self.z[5]*self.z[134]*q5p - rR*self.z[384] - self.z[2]*self.z[385] - self.z[15]*self.z[386] - self.z[128]*self.z[315] -self.z[147]*self.z[351] - self.z[181]*self.z[387] - d1*self.z[6]*self.z[384] - d3*self.z[6]*self.z[375]
        self.z[367] = (self.z[25]*self.z[83]+self.z[92]*self.z[142]+self.z[96]*self.z[145])*self.z[314] + self.z[35]*(self.z[25]*self.z[356]+self.z[83]*self.z[316]+self.z[92]*self.z[339]+self.z[96]*self.z[341]+self.z[142]*self.z[365]+self.z[145]*self.z[366])
        self.z[368] = self.z[34]*self.z[366] + rF*self.z[96]*self.z[343]
        self.z[370] = (self.z[24]*self.z[83]+self.z[92]*self.z[141]+self.z[96]*self.z[144])*self.z[314] + self.z[35]*(self.z[24]*self.z[356]+self.z[83]*self.z[333]+self.z[92]*self.z[346]+self.z[96]*self.z[347]+self.z[141]*self.z[365]+self.z[144]*self.z[366])
        self.z[372] = self.z[34]*self.z[365] + rF*self.z[92]*self.z[343]
        self.z[373] = self.z[1]*self.z[367] + self.z[126]*self.z[317] + self.z[163]*self.z[345] + self.z[178]*self.z[368] + d2*self.z[13]*self.z[356] - self.z[1]*self.z[152]*q3p - self.z[2]*self.z[157]*q3p - self.z[5]*self.z[112]*q5p - self.z[5]*self.z[133]*q5p - rR*self.z[369] - self.z[2]*self.z[370] - self.z[15]*self.z[371] - self.z[130]*self.z[315] -self.z[149]*self.z[351] - self.z[181]*self.z[372] - d1*self.z[6]*self.z[369] - d3*self.z[6]*self.z[356]
        self.z[342] = (self.z[25]*self.z[82]+self.z[93]*self.z[142]+self.z[95]*self.z[145])*self.z[314] + self.z[35]*(self.z[25]*self.z[325]+self.z[82]*self.z[316]+self.z[93]*self.z[339]+self.z[95]*self.z[341]+self.z[142]*self.z[338]+self.z[145]*self.z[340])
        self.z[344] = self.z[34]*self.z[340] + rF*self.z[95]*self.z[343]
        self.z[348] = (self.z[24]*self.z[82]+self.z[93]*self.z[141]+self.z[95]*self.z[144])*self.z[314] + self.z[35]*(self.z[24]*self.z[325]+self.z[82]*self.z[333]+self.z[93]*self.z[346]+self.z[95]*self.z[347]+self.z[141]*self.z[338]+self.z[144]*self.z[340])
        self.z[350] = self.z[34]*self.z[338] + rF*self.z[93]*self.z[343]
        self.z[352] = self.z[1]*self.z[342] + self.z[125]*self.z[317] + self.z[162]*self.z[345] + self.z[178]*self.z[344] + d2*self.z[13]*self.z[325] - self.z[1]*self.z[151]*q3p - self.z[2]*self.z[156]*q3p - self.z[5]*self.z[114]*q5p - self.z[5]*self.z[132]*q5p - self.z[2]*self.z[348] - self.z[15]*self.z[349] - self.z[129]*self.z[315] - self.z[150]*self.z[351] -self.z[181]*self.z[350] - d3*self.z[6]*self.z[325]
        self.z[395] = (self.z[91]*self.z[142]+self.z[97]*self.z[145])*self.z[314] + self.z[35]*(self.z[91]*self.z[339]+self.z[97]*self.z[341]+self.z[142]*self.z[393]+self.z[145]*self.z[394])
        self.z[396] = self.z[34]*self.z[394] + rF*self.z[97]*self.z[343]
        self.z[397] = (self.z[91]*self.z[141]+self.z[97]*self.z[144])*self.z[314] + self.z[35]*(self.z[91]*self.z[346]+self.z[97]*self.z[347]+self.z[141]*self.z[393]+self.z[144]*self.z[394])
        self.z[399] = self.z[34]*self.z[393] + rF*self.z[91]*self.z[343]
        self.z[400] = self.z[1]*self.z[395] + self.z[164]*self.z[345] + self.z[178]*self.z[396] - self.z[1]*self.z[155]*q3p -self.z[2]*self.z[160]*q3p - self.z[2]*self.z[397] - self.z[15]*self.z[398] - self.z[131]*self.z[315] - self.z[148]*self.z[351] - self.z[181]*self.z[399]
        self.z[401] = self.z[89]*(self.z[25]*self.z[314]+self.z[35]*self.z[316])
        self.z[402] = self.z[89]*(self.z[24]*self.z[314]+self.z[35]*self.z[333])
        self.z[403] = self.z[1]*self.z[401] - self.z[1]*self.z[154]*q3p - self.z[2]*self.z[159]*q3p - self.z[2]*self.z[402]
        self.z[404] = self.z[1]*u2*q3p + u3*self.z[388] + u4*self.z[373] + u5*self.z[352] + u7*self.z[400] + u8*self.z[403] - self.z[2]*u1*q3p
        self.z[287] = (self.z[1]*u2-self.z[2]*u1)*q3p
        self.z[407] = self.z[3]*self.z[9]*q7p + self.z[9]*self.z[406] - self.z[4]*self.z[10]*q4p - self.z[10]*self.z[173]*q7p
        self.z[408] = -self.z[3]*self.z[10]*q7p - self.z[4]*self.z[9]*q4p - self.z[9]*self.z[173]*q7p - self.z[10]*self.z[406]
        self.z[409] = self.z[11]*self.z[76]*q8p + self.z[11]*self.z[407] + self.z[12]*self.z[355] - self.z[12]*self.z[174]*q8p
        self.z[410] = self.z[11]*self.z[174]*q8p + self.z[12]*self.z[76]*q8p + self.z[12]*self.z[407] - self.z[11]*self.z[355]
        self.z[411] = self.z[2]*self.z[153]*q3p + self.z[113]*self.z[355] + self.z[127]*self.z[407] + self.z[128]*self.z[408] +self.z[134]*self.z[355] + self.z[161]*self.z[409] + self.z[176]*self.z[386] + self.z[179]*self.z[383] + rR*self.z[3]*self.z[294] + d1*self.z[3]*self.z[405] + d1*self.z[76]*self.z[384] + d2*self.z[174]*self.z[375] + d3*self.z[76]*self.z[375] - self.z[1]*self.z[158]*q3p - self.z[4]*self.z[100]*q4p - self.z[4]*self.z[110]*q4p - self.z[1]*self.z[385] -self.z[2]*self.z[382] - self.z[147]*self.z[410] - self.z[182]*self.z[387]
        self.z[413] = self.z[2]*self.z[152]*q3p + self.z[112]*self.z[355] + self.z[126]*self.z[407] + self.z[130]*self.z[408] +self.z[133]*self.z[355] + self.z[163]*self.z[409] + self.z[176]*self.z[371] + self.z[179]*self.z[368] + rR*self.z[3]*self.z[304] + d1*self.z[3]*self.z[412] + d1*self.z[76]*self.z[369] + d2*self.z[174]*self.z[356] + d3*self.z[76]*self.z[356] - self.z[1]*self.z[157]*q3p - self.z[4]*self.z[101]*q4p - self.z[4]*self.z[111]*q4p - self.z[1]*self.z[370] -self.z[2]*self.z[367] - self.z[149]*self.z[410] - self.z[182]*self.z[372]
        self.z[414] = self.z[2]*self.z[151]*q3p + self.z[114]*self.z[355] + self.z[125]*self.z[407] + self.z[129]*self.z[408] +self.z[132]*self.z[355] + self.z[162]*self.z[409] + self.z[176]*self.z[349] + self.z[179]*self.z[344] + d2*self.z[174]*self.z[325] + d3*self.z[76]*self.z[325] - self.z[1]*self.z[156]*q3p - self.z[1]*self.z[348] - self.z[2]*self.z[342] -self.z[150]*self.z[410] - self.z[182]*self.z[350]
        self.z[415] = self.z[2]*self.z[155]*q3p + self.z[131]*self.z[408] + self.z[164]*self.z[409] + self.z[176]*self.z[398] +self.z[179]*self.z[396] - self.z[1]*self.z[160]*q3p - self.z[1]*self.z[397] - self.z[2]*self.z[395] - self.z[148]*self.z[410] -self.z[182]*self.z[399]
        self.z[416] = self.z[2]*self.z[154]*q3p - self.z[1]*self.z[159]*q3p - self.z[1]*self.z[402] - self.z[2]*self.z[401]
        self.z[417] = u3*self.z[411] + u4*self.z[413] + u5*self.z[414] + u7*self.z[415] + u8*self.z[416] - self.z[1]*u1*q3p - self.z[2]*u2*q3p
        self.z[436] = (self.z[262]*self.z[310]+self.z[219]*self.z[236]*self.z[404]-self.z[250]*self.z[287]-self.z[270]*self.z[433]-self.z[219]*self.z[265]*self.z[417])/self.z[234]
        self.z[803] = self.z[744] + self.z[747]*self.z[275] + self.z[281]*self.z[748] + self.z[283]*self.z[749]
        self.z[434] = (self.z[257]*self.z[310]+self.z[228]*self.z[236]*self.z[404]-self.z[243]*self.z[287]-self.z[268]*self.z[433]-self.z[228]*self.z[265]*self.z[417])/self.z[234]
        self.z[806] = self.z[281]*self.z[759] + self.z[283]*self.z[766] - self.z[286]*self.z[771]
        self.z[438] = (self.z[263]*self.z[310]+self.z[267]*self.z[417]-self.z[251]*self.z[287]-self.z[264]*self.z[404]-self.z[272]*self.z[433])/self.z[234]
        self.z[804] = self.z[745] + self.z[281]*self.z[751] + self.z[283]*self.z[752] - self.z[747]*self.z[278]
        self.z[435] = (self.z[261]*self.z[310]+self.z[223]*self.z[236]*self.z[404]-self.z[248]*self.z[287]-self.z[269]*self.z[433]-self.z[223]*self.z[265]*self.z[417])/self.z[234]
        self.z[811] = self.z[807] + self.z[799]*self.z[437] + self.z[801]*self.z[436] + self.z[803]*self.z[434] + self.z[806]*self.z[438] - self.z[804]*self.z[435]
        self.z[611] = ic12*self.z[44]
        self.z[620] = ic22*self.z[39]
        self.z[632] = ic23*self.z[46]
        self.z[646] = id22*self.z[39]
        self.z[608] = ic11*self.z[44]
        self.z[616] = ic12*self.z[39]
        self.z[629] = ic31*self.z[46]
        self.z[614] = ic31*self.z[44]
        self.z[624] = ic23*self.z[39]
        self.z[635] = ic33*self.z[46]
        self.z[644] = id11*self.z[52]
        self.z[652] = id33*self.z[54]
        self.z[662] = ie11*self.z[80]
        self.z[676] = ie12*self.z[83]
        self.z[689] = ie31*self.z[87]
        self.z[667] = ie12*self.z[80]
        self.z[680] = ie22*self.z[83]
        self.z[694] = ie23*self.z[87]
        self.z[714] = if22*self.z[83]
        self.z[672] = ie31*self.z[80]
        self.z[684] = ie23*self.z[83]
        self.z[699] = ie33*self.z[87]
        self.z[710] = if11*self.z[92]
        self.z[720] = if33*self.z[96]
        self.z[727] = self.z[39]*self.z[611] + self.z[39]*self.z[620] + self.z[39]*self.z[632] + self.z[39]*self.z[646] + self.z[44]*self.z[608] + self.z[44]*self.z[616] + self.z[44]*self.z[629] + self.z[46]*self.z[614] + self.z[46]*self.z[624] + self.z[46]*self.z[635] + self.z[52]*self.z[644] + self.z[54]*self.z[652] + self.z[80]*self.z[662] + self.z[80]*self.z[676] + self.z[80]*self.z[689] + self.z[83]*self.z[667] + self.z[83]*self.z[680] + self.z[83]*self.z[694] + self.z[83]*self.z[714] + self.z[87]*self.z[672] + self.z[87]*self.z[684] + self.z[87]*self.z[699] + self.z[92]*self.z[710] + self.z[96]*self.z[720] + massd*(pow(self.z[98],2)+pow(self.z[101],2)) - massc*(2*self.z[5]*self.z[98]*self.z[102]-2*self.z[101]*self.z[106]-pow(self.z[98],2)-pow(self.z[101],2)-pow(self.z[102],2)-pow(self.z[106],2)-pow(self.z[107],2)-2*self.z[6]*self.z[98]*self.z[107]) - masse*(2*self.z[13]*self.z[98]*self.z[116]-2*self.z[101]*self.z[111]-2*self.z[112]*self.z[123]-pow(self.z[98],2)-pow(self.z[101],2)-pow(self.z[111],2)-pow(self.z[112],2)-pow(self.z[116],2)-pow(self.z[120],2)-pow(self.z[123],2)-2*self.z[6]*self.z[98]*self.z[112]-2*self.z[6]*self.z[98]*self.z[123]-2*self.z[9]*self.z[101]*self.z[120]-2*self.z[9]*self.z[111]*self.z[120]-2*self.z[10]*self.z[101]*self.z[116]-2*self.z[10]*self.z[111]*self.z[116]-2*self.z[15]*self.z[98]*self.z[120]) - massf*(2*self.z[13]*self.z[98]*self.z[126]-2*self.z[101]*self.z[111]-2*self.z[112]*self.z[133]-pow(self.z[98],2)-pow(self.z[101],2)-pow(self.z[111],2)-pow(self.z[112],2)-pow(self.z[126],2)-pow(self.z[130],2)-pow(self.z[133],2)-2*self.z[6]*self.z[98]*self.z[112]-2*self.z[6]*self.z[98]*self.z[133]-2*self.z[9]*self.z[101]*self.z[130]-2*self.z[9]*self.z[111]*self.z[130]-2*self.z[10]*self.z[101]*self.z[126]-2*self.z[10]*self.z[111]*self.z[126]-2*self.z[15]*self.z[98]*self.z[130])
        self.z[755] = self.z[40]*self.z[611] + self.z[40]*self.z[620] + self.z[40]*self.z[632] + self.z[40]*self.z[646] + self.z[43]*self.z[608] + self.z[43]*self.z[616] + self.z[43]*self.z[629] + self.z[45]*self.z[614] + self.z[45]*self.z[624] + self.z[45]*self.z[635] + self.z[51]*self.z[644] + self.z[53]*self.z[652] + self.z[78]*self.z[662] + self.z[78]*self.z[676] + self.z[78]*self.z[689] + self.z[84]*self.z[667] + self.z[84]*self.z[680] + self.z[84]*self.z[694] + self.z[84]*self.z[714] + self.z[85]*self.z[672] + self.z[85]*self.z[684] + self.z[85]*self.z[699] + self.z[90]*self.z[710] + self.z[94]*self.z[720] + massd*(self.z[98]*self.z[99]+self.z[100]*self.z[101]) - massc*(self.z[5]*self.z[98]*self.z[103]+self.z[5]*self.z[99]*self.z[102]-self.z[98]*self.z[99]-self.z[100]*self.z[101]-self.z[100]*self.z[106]-self.z[101]*self.z[105]-self.z[102]*self.z[103]-self.z[105]*self.z[106]-self.z[107]*self.z[108]-self.z[6]*self.z[98]*self.z[108]-self.z[6]*self.z[99]*self.z[107]) - masse*(self.z[13]*self.z[98]*self.z[117]+self.z[13]*self.z[99]*self.z[116]-self.z[98]*self.z[99]-self.z[100]*self.z[101]-self.z[100]*self.z[111]-self.z[101]*self.z[110]-self.z[110]*self.z[111]-self.z[112]*self.z[113]-self.z[112]*self.z[124]-self.z[113]*self.z[123]-self.z[116]*self.z[117]-self.z[118]*self.z[120]-self.z[123]*self.z[124]-self.z[6]*self.z[98]*self.z[113]-self.z[6]*self.z[98]*self.z[124]-self.z[6]*self.z[99]*self.z[112]-self.z[6]*self.z[99]*self.z[123]-self.z[9]*self.z[100]*self.z[120]-self.z[9]*self.z[101]*self.z[118]-self.z[9]*self.z[110]*self.z[120]-self.z[9]*self.z[111]*self.z[118]-self.z[10]*self.z[100]*self.z[116]-self.z[10]*self.z[101]*self.z[117]-self.z[10]*self.z[110]*self.z[116]-self.z[10]*self.z[111]*self.z[117]-self.z[15]*self.z[98]*self.z[118]-self.z[15]*self.z[99]*self.z[120]) - massf*(self.z[13]*self.z[98]*self.z[127]+self.z[13]*self.z[99]*self.z[126]-self.z[98]*self.z[99]-self.z[100]*self.z[101]-self.z[100]*self.z[111]-self.z[101]*self.z[110]-self.z[110]*self.z[111]-self.z[112]*self.z[113]-self.z[112]*self.z[134]-self.z[113]*self.z[133]-self.z[126]*self.z[127]-self.z[128]*self.z[130]-self.z[133]*self.z[134]-self.z[6]*self.z[98]*self.z[113]-self.z[6]*self.z[98]*self.z[134]-self.z[6]*self.z[99]*self.z[112]-self.z[6]*self.z[99]*self.z[133]-self.z[9]*self.z[100]*self.z[130]-self.z[9]*self.z[101]*self.z[128]-self.z[9]*self.z[110]*self.z[130]-self.z[9]*self.z[111]*self.z[128]-self.z[10]*self.z[100]*self.z[126]-self.z[10]*self.z[101]*self.z[127]-self.z[10]*self.z[110]*self.z[126]-self.z[10]*self.z[111]*self.z[127]-self.z[15]*self.z[98]*self.z[128]-self.z[15]*self.z[99]*self.z[130])
        self.z[762] = self.z[42]*self.z[646] + self.z[81]*self.z[662] + self.z[81]*self.z[676] + self.z[81]*self.z[689] + self.z[82]*self.z[667] + self.z[82]*self.z[680] + self.z[82]*self.z[694] + self.z[82]*self.z[714] + self.z[86]*self.z[672] + self.z[86]*self.z[684] + self.z[86]*self.z[699] + self.z[93]*self.z[710] + self.z[95]*self.z[720] + self.z[42]*(self.z[611]+self.z[620]+self.z[632]) - massc*(self.z[5]*self.z[98]*self.z[104]-self.z[102]*self.z[104]-self.z[107]*self.z[109]-self.z[6]*self.z[98]*self.z[109]) - masse*(self.z[13]*self.z[98]*self.z[115]-self.z[112]*self.z[114]-self.z[112]*self.z[122]-self.z[114]*self.z[123]-self.z[115]*self.z[116]-self.z[119]*self.z[120]-self.z[122]*self.z[123]-self.z[6]*self.z[98]*self.z[114]-self.z[6]*self.z[98]*self.z[122]-self.z[9]*self.z[101]*self.z[119]-self.z[9]*self.z[111]*self.z[119]-self.z[10]*self.z[101]*self.z[115]-self.z[10]*self.z[111]*self.z[115]-self.z[15]*self.z[98]*self.z[119]) - massf*(self.z[13]*self.z[98]*self.z[125]-self.z[112]*self.z[114]-self.z[112]*self.z[132]-self.z[114]*self.z[133]-self.z[125]*self.z[126]-self.z[129]*self.z[130]-self.z[132]*self.z[133]-self.z[6]*self.z[98]*self.z[114]-self.z[6]*self.z[98]*self.z[132]-self.z[9]*self.z[101]*self.z[129]-self.z[9]*self.z[111]*self.z[129]-self.z[10]*self.z[101]*self.z[125]-self.z[10]*self.z[111]*self.z[125]-self.z[15]*self.z[98]*self.z[129])
        self.z[769] = self.z[89]*self.z[714]
        self.z[774] = self.z[727] + self.z[274]*self.z[733] + self.z[280]*self.z[755] + self.z[282]*self.z[762] - self.z[277]*self.z[732] - self.z[285]*self.z[769]
        self.z[782] = self.z[774] + self.z[274]*self.z[780] + self.z[280]*self.z[773] + self.z[282]*self.z[775] - self.z[277]*self.z[779] - self.z[285]*self.z[778]
        self.z[735] = self.z[47]*self.z[646]
        self.z[786] = self.z[735] + self.z[276]*self.z[732] - self.z[273]*self.z[733] - self.z[279]*self.z[755] - self.z[284]*self.z[769]
        self.z[795] = self.z[786] + self.z[274]*self.z[790] + self.z[280]*self.z[787] + self.z[282]*self.z[788] - self.z[277]*self.z[792] - self.z[285]*self.z[793]
        self.z[816] = self.z[782]*self.z[797] - self.z[784]*self.z[795]
        self.z[815] = self.z[782]*self.z[796] - self.z[783]*self.z[795]
        self.z[165] = l6*self.z[39]
        self.z[166] = l6*self.z[40]
        self.z[167] = l6*self.z[42]
        self.z[168] = l5*self.z[45] - l6*self.z[43]
        self.z[169] = l5*self.z[46] - l6*self.z[44]
        self.z[170] = l5*self.z[39]
        self.z[171] = l5*self.z[40]
        self.z[172] = l5*self.z[42]
        self.z[519] = self.z[46] + self.z[45]*self.z[280]
        self.z[521] = self.z[40]*self.z[279]
        self.z[522] = self.z[45]*self.z[279]
        self.z[525] = self.z[45]*self.z[281]
        self.z[528] = self.z[47] - self.z[40]*self.z[279]
        self.z[535] = self.z[87] + self.z[85]*self.z[280] + self.z[86]*self.z[282]
        self.z[538] = self.z[85]*self.z[279]
        self.z[540] = self.z[88] + self.z[85]*self.z[281] + self.z[86]*self.z[283]
        self.z[542] = -self.z[98] - self.z[99]*self.z[280]
        self.z[543] = self.z[101] + self.z[100]*self.z[280]
        self.z[544] = self.z[102] + self.z[103]*self.z[280] + self.z[104]*self.z[282]
        self.z[545] = self.z[106] + self.z[105]*self.z[280]
        self.z[546] = -self.z[107] - self.z[108]*self.z[280] - self.z[109]*self.z[282]
        self.z[547] = self.z[99]*self.z[279]
        self.z[548] = self.z[100]*self.z[279]
        self.z[549] = self.z[103]*self.z[279]
        self.z[550] = self.z[105]*self.z[279]
        self.z[551] = self.z[108]*self.z[279]
        self.z[552] = self.z[99]*self.z[281]
        self.z[553] = self.z[100]*self.z[281]
        self.z[554] = self.z[103]*self.z[281] + self.z[104]*self.z[283]
        self.z[555] = self.z[105]*self.z[281]
        self.z[556] = -self.z[108]*self.z[281] - self.z[109]*self.z[283]
        self.z[557] = self.z[111] + self.z[110]*self.z[280]
        self.z[558] = -self.z[112] - self.z[113]*self.z[280] - self.z[114]*self.z[282]
        self.z[559] = self.z[116] + self.z[115]*self.z[282] + self.z[117]*self.z[280]
        self.z[560] = self.z[120] + self.z[118]*self.z[280] + self.z[119]*self.z[282]
        self.z[561] = -self.z[123] - self.z[122]*self.z[282] - self.z[124]*self.z[280]
        self.z[562] = self.z[110]*self.z[279]
        self.z[563] = self.z[113]*self.z[279]
        self.z[564] = self.z[117]*self.z[279]
        self.z[565] = self.z[118]*self.z[279]
        self.z[566] = self.z[124]*self.z[279]
        self.z[567] = self.z[121] + self.z[118]*self.z[281] + self.z[119]*self.z[283]
        self.z[568] = self.z[110]*self.z[281]
        self.z[569] = -self.z[113]*self.z[281] - self.z[114]*self.z[283]
        self.z[570] = self.z[115]*self.z[283] + self.z[117]*self.z[281]
        self.z[571] = -self.z[122]*self.z[283] - self.z[124]*self.z[281]
        self.z[572] = self.z[126] + self.z[125]*self.z[282] + self.z[127]*self.z[280]
        self.z[573] = self.z[130] + self.z[128]*self.z[280] + self.z[129]*self.z[282]
        self.z[574] = -self.z[133] - self.z[132]*self.z[282] - self.z[134]*self.z[280]
        self.z[575] = self.z[127]*self.z[279]
        self.z[576] = self.z[128]*self.z[279]
        self.z[577] = self.z[134]*self.z[279]
        self.z[578] = self.z[131] + self.z[128]*self.z[281] + self.z[129]*self.z[283]
        self.z[579] = self.z[125]*self.z[283] + self.z[127]*self.z[281]
        self.z[580] = -self.z[132]*self.z[283] - self.z[134]*self.z[281]
        self.z[584] = self.z[38] + self.z[37]*self.z[280]
        self.z[587] = self.z[37]*self.z[279]
        self.z[589] = self.z[37]*self.z[281]
        self.z[592] = self.z[165] + self.z[166]*self.z[280] + self.z[167]*self.z[282]
        self.z[593] = self.z[169] + self.z[168]*self.z[280]
        self.z[594] = -self.z[170] - self.z[171]*self.z[280] - self.z[172]*self.z[282]
        self.z[595] = self.z[166]*self.z[279]
        self.z[596] = self.z[168]*self.z[279]
        self.z[597] = self.z[171]*self.z[279]
        self.z[598] = self.z[166]*self.z[281] + self.z[167]*self.z[283]
        self.z[599] = self.z[168]*self.z[281]
        self.z[600] = -self.z[171]*self.z[281] - self.z[172]*self.z[283]
        self.z[734] = self.z[39]*self.z[612] + self.z[39]*self.z[623] + self.z[39]*self.z[633] + self.z[39]*self.z[638] + self.z[39]*self.z[650] + self.z[39]*self.z[655] + self.z[44]*self.z[609] + self.z[44]*self.z[619] + self.z[44]*self.z[630] + self.z[44]*self.z[639] + self.z[46]*self.z[615] + self.z[46]*self.z[627] + self.z[46]*self.z[636] + self.z[46]*self.z[637] + self.z[52]*self.z[645] + self.z[52]*self.z[656] + self.z[54]*self.z[653] + self.z[54]*self.z[654] + self.z[80]*self.z[664] + self.z[80]*self.z[678] + self.z[80]*self.z[691] + self.z[80]*self.z[704] + self.z[83]*self.z[669] + self.z[83]*self.z[682] + self.z[83]*self.z[696] + self.z[83]*self.z[703] + self.z[83]*self.z[717] + self.z[83]*self.z[724] + self.z[87]*self.z[674] + self.z[87]*self.z[686] + self.z[87]*self.z[701] + self.z[87]*self.z[702] + self.z[92]*self.z[712] + self.z[92]*self.z[725] + self.z[96]*self.z[722] + self.z[96]*self.z[723] + massc*(self.z[101]*self.z[479]+self.z[101]*self.z[482]+self.z[102]*self.z[481]+self.z[106]*self.z[479]+self.z[106]*self.z[482]+self.z[5]*self.z[102]*self.z[478]-self.z[98]*self.z[478]-self.z[107]*self.z[483]-self.z[5]*self.z[98]*self.z[481]-self.z[5]*self.z[107]*self.z[480]-self.z[6]*self.z[98]*self.z[483]-self.z[6]*self.z[102]*self.z[480]-self.z[6]*self.z[107]*self.z[478]) - massd*(self.z[98]*self.z[478]-self.z[101]*self.z[479]) - masse*(self.z[98]*self.z[478]+self.z[112]*self.z[500]+self.z[112]*self.z[503]+self.z[123]*self.z[500]+self.z[123]*self.z[503]+self.z[5]*self.z[98]*self.z[498]+self.z[5]*self.z[112]*self.z[480]+self.z[5]*self.z[123]*self.z[480]+self.z[6]*self.z[98]*self.z[500]+self.z[6]*self.z[98]*self.z[503]+self.z[6]*self.z[112]*self.z[478]+self.z[6]*self.z[123]*self.z[478]+self.z[10]*self.z[120]*self.z[498]+self.z[13]*self.z[98]*self.z[501]+self.z[14]*self.z[116]*self.z[480]+self.z[15]*self.z[120]*self.z[478]-self.z[101]*self.z[479]-self.z[101]*self.z[499]-self.z[111]*self.z[479]-self.z[111]*self.z[499]-self.z[116]*self.z[501]-self.z[120]*self.z[502]-self.z[9]*self.z[101]*self.z[502]-self.z[9]*self.z[111]*self.z[502]-self.z[9]*self.z[116]*self.z[498]-self.z[9]*self.z[120]*self.z[479]-self.z[9]*self.z[120]*self.z[499]-self.z[10]*self.z[101]*self.z[501]-self.z[10]*self.z[111]*self.z[501]-self.z[10]*self.z[116]*self.z[479]-self.z[10]*self.z[116]*self.z[499]-self.z[13]*self.z[116]*self.z[478]-self.z[15]*self.z[98]*self.z[502]-self.z[16]*self.z[120]*self.z[480]) -massf*(self.z[98]*self.z[478]+self.z[112]*self.z[500]+self.z[112]*self.z[512]+self.z[133]*self.z[500]+self.z[133]*self.z[512]+self.z[5]*self.z[98]*self.z[498]+self.z[5]*self.z[112]*self.z[480]+self.z[5]*self.z[133]*self.z[480]+self.z[6]*self.z[98]*self.z[500]+self.z[6]*self.z[98]*self.z[512]+self.z[6]*self.z[112]*self.z[478]+self.z[6]*self.z[133]*self.z[478]+self.z[10]*self.z[130]*self.z[498]+self.z[13]*self.z[98]*self.z[510]+self.z[14]*self.z[126]*self.z[480]+self.z[15]*self.z[130]*self.z[478]-self.z[101]*self.z[479]-self.z[101]*self.z[499]-self.z[111]*self.z[479]-self.z[111]*self.z[499]-self.z[126]*self.z[510]-self.z[130]*self.z[511]-self.z[9]*self.z[101]*self.z[511]-self.z[9]*self.z[111]*self.z[511]-self.z[9]*self.z[126]*self.z[498]-self.z[9]*self.z[130]*self.z[479]-self.z[9]*self.z[130]*self.z[499]-self.z[10]*self.z[101]*self.z[510]-self.z[10]*self.z[111]*self.z[510]-self.z[10]*self.z[126]*self.z[479]-self.z[10]*self.z[126]*self.z[499]-self.z[13]*self.z[126]*self.z[478]-self.z[15]*self.z[98]*self.z[511]-self.z[16]*self.z[130]*self.z[480])
        self.z[739] = self.z[47]*(self.z[650]+self.z[655])
        self.z[741] = self.z[79]*self.z[662] + self.z[79]*self.z[676] + self.z[79]*self.z[689] + self.z[88]*self.z[672] + self.z[88]*self.z[684] + self.z[88]*self.z[699] + self.z[91]*self.z[710] + self.z[97]*self.z[720] + masse*self.z[121]*(self.z[120]+self.z[9]*self.z[101]+self.z[9]*self.z[111]+self.z[15]*self.z[98]) + massf*self.z[131]*(self.z[130]+self.z[9]*self.z[101]+self.z[9]*self.z[111]+self.z[15]*self.z[98])
        self.z[743] = self.z[79]*self.z[661] + self.z[79]*self.z[690] + self.z[88]*self.z[671] + self.z[88]*self.z[700] + self.z[91]*self.z[709] + self.z[97]*self.z[721] + masse*pow(self.z[121],2) + massf*pow(self.z[131],2)
        self.z[781] = self.z[734] + self.z[274]*self.z[750] + self.z[280]*self.z[760] + self.z[282]*self.z[767] - self.z[277]*self.z[753] - self.z[285]*self.z[772]
        self.z[785] = self.z[781] + self.z[773]*self.z[436] + self.z[775]*self.z[437] + self.z[778]*self.z[438] + self.z[780]*self.z[434] - self.z[779]*self.z[435]
        self.z[794] = self.z[739] + self.z[276]*self.z[753] - self.z[273]*self.z[750] - self.z[279]*self.z[760] - self.z[284]*self.z[772]
        self.z[798] = self.z[794] + self.z[787]*self.z[436] + self.z[788]*self.z[437] + self.z[790]*self.z[434] + self.z[793]*self.z[438] - self.z[792]*self.z[435]
        self.z[800] = self.z[741] + self.z[275]*self.z[733] + self.z[281]*self.z[755] + self.z[283]*self.z[762] - self.z[278]*self.z[732] - self.z[286]*self.z[769]
        self.z[802] = self.z[743] + self.z[275]*self.z[744] + self.z[281]*self.z[758] + self.z[283]*self.z[765] - self.z[278]*self.z[745]
        self.z[805] = self.z[281]*self.z[757] + self.z[283]*self.z[764]
        self.z[808] = self.z[800] + self.z[274]*self.z[803] + self.z[280]*self.z[801] + self.z[282]*self.z[799] - self.z[277]*self.z[804] - self.z[285]*self.z[806]
        self.z[809] = self.z[805] + self.z[276]*self.z[804] - self.z[273]*self.z[803] - self.z[279]*self.z[801] - self.z[284]*self.z[806]
        self.z[810] = self.z[802] + self.z[275]*self.z[803] + self.z[281]*self.z[801] + self.z[283]*self.z[799] - self.z[278]*self.z[804] - self.z[286]*self.z[806]
        self.z[818] = self.z[809]*self.z[816] - self.z[808]*self.z[817] - self.z[810]*self.z[815]
        self.z[819] = self.z[796]*self.z[810] - self.z[797]*self.z[809]
        self.z[820] = self.z[795]*self.z[810] - self.z[797]*self.z[808]
        self.z[821] = self.z[795]*self.z[809] - self.z[796]*self.z[808]
        self.z[822] = self.z[783]*self.z[810] - self.z[784]*self.z[809]
        self.z[823] = self.z[782]*self.z[810] - self.z[784]*self.z[808]
        self.z[824] = self.z[782]*self.z[809] - self.z[783]*self.z[808]
        T4 = Tphi
        T6 = TthetaR
        T7 = Tdelta
        self.z[603] = T4*self.z[589] + T7*self.z[540] + self.z[514]*self.z[4]*self.z[553] + self.z[513]*(self.z[4]*self.z[553]+self.z[4]*self.z[555]+self.z[29]*self.z[556]-self.z[77]*self.z[554]) + self.z[515]*(self.z[4]*self.z[553]+self.z[4]*self.z[568]+self.z[23]*self.z[570]+self.z[26]*self.z[567]+self.z[29]*self.z[569]+self.z[29]*self.z[571]) + self.z[516]*(self.z[4]*self.z[553]+self.z[4]*self.z[568]+self.z[23]*self.z[579]+self.z[26]*self.z[578]+self.z[29]*self.z[569]+self.z[29]*self.z[580]) - T7*self.z[525] - Fphi*(self.z[278]+self.z[2]*self.z[552]-self.z[18]*self.z[553]-self.z[18]*self.z[599]-self.z[28]*self.z[600]-self.z[74]*self.z[598])
        self.z[814] = self.z[811] - self.z[603]
        self.z[601] = T4*self.z[584] + T7*self.z[535] + self.z[514]*self.z[4]*self.z[543] + self.z[513]*(self.z[4]*self.z[543]+self.z[4]*self.z[545]+self.z[29]*self.z[546]-self.z[77]*self.z[544]) + self.z[515]*(self.z[4]*self.z[543]+self.z[4]*self.z[557]+self.z[23]*self.z[559]+self.z[26]*self.z[560]+self.z[29]*self.z[558]+self.z[29]*self.z[561]) + self.z[516]*(self.z[4]*self.z[543]+self.z[4]*self.z[557]+self.z[23]*self.z[572]+self.z[26]*self.z[573]+self.z[29]*self.z[558]+self.z[29]*self.z[574]) - T7*self.z[519] - Fphi*(self.z[277]-self.z[2]*self.z[542]-self.z[18]*self.z[543]-self.z[18]*self.z[593]-self.z[28]*self.z[594]-self.z[74]*self.z[592])
        self.z[812] = self.z[785] - self.z[601]
        self.z[602] = T6*self.z[521] + T6*self.z[528] + T7*self.z[522] + Fphi*(self.z[276]+self.z[2]*self.z[547]+self.z[28]*self.z[597]-self.z[18]*self.z[548]-self.z[18]*self.z[596]-self.z[74]*self.z[595]) - T4*self.z[587] - T7*self.z[538] -self.z[514]*self.z[4]*self.z[548] - self.z[513]*(self.z[4]*self.z[548]+self.z[4]*self.z[550]-self.z[29]*self.z[551]-self.z[77]*self.z[549]) - self.z[515]*(self.z[4]*self.z[548]+self.z[4]*self.z[562]+self.z[23]*self.z[564]+self.z[26]*self.z[565]-self.z[29]*self.z[563]-self.z[29]*self.z[566]) - self.z[516]*(self.z[4]*self.z[548]+self.z[4]*self.z[562]+self.z[23]*self.z[575]+self.z[26]*self.z[576]-self.z[29]*self.z[563]-self.z[29]*self.z[577])
        self.z[813] = self.z[798] - self.z[602]
        self.z[825] = (self.z[817]*self.z[814]+self.z[819]*self.z[812]-self.z[822]*self.z[813])/self.z[818]
        u4p = self.z[825]
        self.z[826] = (self.z[816]*self.z[814]+self.z[820]*self.z[812]-self.z[823]*self.z[813])/self.z[818]
        u6p = -self.z[826]
        self.z[827] = (self.z[815]*self.z[814]+self.z[821]*self.z[812]-self.z[824]*self.z[813])/self.z[818]
        u7p = self.z[827]

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
        u[0] = 0
        u[1] = 0
        u[2] = 0
        u[3] = 0
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
        u1 = self.z[274]*u4 + self.z[275]*u7 - self.z[273]*u6
        u2 = self.z[276]*u6 - self.z[277]*u4 - self.z[278]*u7
        u3 = self.z[280]*u4 + self.z[281]*u7 - self.z[279]*u6
        u5 = self.z[282]*u4 + self.z[283]*u7
        u8 = -self.z[284]*u6 - self.z[285]*u4 - self.z[286]*u7

        # calculate the outputs
        self.z[518] = self.z[39] + self.z[40]*self.z[280] + self.z[42]*self.z[282]
        self.z[524] = self.z[40]*self.z[281] + self.z[42]*self.z[283]
        self.z[533] = self.z[80] + self.z[78]*self.z[280] + self.z[81]*self.z[282]
        self.z[534] = self.z[83] + self.z[82]*self.z[282] + self.z[84]*self.z[280]
        self.z[539] = self.z[79] + self.z[78]*self.z[281] + self.z[81]*self.z[283]
        self.z[541] = self.z[82]*self.z[283] + self.z[84]*self.z[281]
        q9 = q1 + d1*self.z[71] + d2*self.z[27] + d3*self.z[21] - rR*self.z[19] - self.z[24]*self.z[34]
        q10 = q2 + rR*self.z[20] + d1*self.z[74] + d2*self.z[28] + d3*self.z[22] - self.z[25]*self.z[34]

        # plug in the derivatives for returning
        y = zeros(len(self.outputNames))
        for i, name in enumerate(self.outputNames):
            exec('y[' + str(i) + '] = ' + name)

        return y
