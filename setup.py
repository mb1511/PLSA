from setuptools import setup

def readme():
    with open('README.md') as s:
        return s.read()

setup(
    name='plsa',
    version='0.1',
    description='Pan-Locus Proteome Comparison',
    long_description=readme(),
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
    ],
    url='https://github.com/mb1511/PLSA',
    author='Matt Brewer',
    author_email='mb1511@bristol.ac.uk',
    license='GPLv3',
    packages=['plsa'],
    install_requires=['dill', 'ete3', 'numpy'],
    zip_safe=False)

