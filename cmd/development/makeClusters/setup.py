from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='makeClusters',
        python_requires='>=3',
        version='210105',
        description='Building sequence clusters and core/flexible genomes',
        url='https://github.com/apsteinberg/mcorr',
        license='MIT',
        author='Asher Preska Steinberg',
        author_email='apsteinberg@nyu.edu',
        packages=['seqClusters'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['makeClusters=seqClusters.cli:main'],
            }
      )
