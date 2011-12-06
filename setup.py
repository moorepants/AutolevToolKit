from distutils.core import setup

setup(
    name='AutolevToolKit',
    version='0.1.0dev',
    author=['Jason Keith Moore', 'Dale Lukas Peterson'],
    author_email=['moorepants@gmail.com', 'hazelnusse@gmail.com'],
    packages=['altk', 'altk.test', 'models'],
    url='http://github.com/moorepants/AutolevToolKit',
    license='LICENSE.txt',
    description='Parses Autolev 4.1 input and output files for use in other languages.',
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
