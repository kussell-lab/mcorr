from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='mcorr',
        python_requires='>=3',
        version='20180506',
        description='Inferring recombination rates from correlation profiles',
        url='https://github.com/apsteinberg/mcorr',
        author='Mingzhi Lin, Asher Preska Steinberg',
        author_email='mingzhi9@gmail.com, apsteinberg@nyu.edu',
        license='MIT',
        packages=['mcorr'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['mcorr-fit=mcorr.cli:main', 'mcorrFitOne=mcorr.singleFit:main',
                                 'mcorrFitCompare=mcorr.FitComparison:main'],
            },
        zip_safe=False)
