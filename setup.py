from distutils.core import setup

setup(
    name='alparse',
    version='0.1.0dev',
    author=['Dale Lukas Peterson', 'Jason Keith Moore'],
    author_email=['hazelnusse@gmail.com', 'moorepants@gmail.com'],
    packages=['alparse', 'alparse.test', 'models'],
    url='http://github.com/hazelnusse/alparse',
    license='LICENSE.txt',
    description='Parses Autolev 4.1 output for use in other languages.',
    long_description=open('README.rst').read(),
    classifiers=['Programming Language :: Python',
                 'Programming Language :: Python :: 2.7',
                 'Operating System :: OS Independent',
                 'Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: BSD License',
                 'Natural Language :: English',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Physics']
)
