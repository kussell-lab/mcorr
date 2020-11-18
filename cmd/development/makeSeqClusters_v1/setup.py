from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='makeSeqClusters',
        python_requires='>=3',
        version='201117',
        description='Building sequence clusters and core/flexible genomes',
        url='https://github.com/apsteinberg/mcorr',
        author='Asher Preska Steinberg',
        author_email='apsteinberg@nyu.edu',
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['makeSeqClusters=makeSeqClusters:main'],
            }
      )
