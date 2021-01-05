from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='clusterSequences',
        python_requires='>=3',
        version='210105',
        description='cluster sequences based on pairwise distances from mcorr-pair/mcorr-pair-sync',
        url='https://github.com/apsteinberg/mcorr',
        license='MIT',
        author='Asher Preska Steinberg',
        author_email='apsteinberg@nyu.edu',
        packages=['bin'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['clusterSequences=bin.clusterSequences:main'],
            }
      )
