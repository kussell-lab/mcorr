from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='mcorr',
        python_requires='>=3',
        version='20180314',
        description='Inferring recombination rates from correlation profiles',
        url='https://github.com/kussell-lab/mcorr',
        author='Mingzhi Lin',
        author_email='mingzhi9@gmail.com',
        license='MIT',
        packages=['mcorr'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['mcorr-fit=mcorr.cli:main'],
            },
        zip_safe=False)
