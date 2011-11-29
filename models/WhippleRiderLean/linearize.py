import sys
sys.path.append('/media/Data/Documents/School/Autolev/autolev_toolkit')
import autolev_toolkit as altk

states = ['q' + str(x + 1) for x in range(9)] + ['u4', 'u6', 'u7', 'u9']

inputs = ['T4', 'T6', 'T7', 'T9']

outputs = (['q' + str(x + 1) for x in range(9)] +
           ['u' + str(x + 1) for x in range(9)])

altk.write_linearization(('A', 'B', 'C', 'D'), states, inputs, outputs,
        filename='linearize.al')
