from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='getDivergence',
        python_requires='>=3',
        version='201117',
        description='Collect results for many sequence clusters from mcorr-fit',
        url='https://github.com/apsteinberg/mcorr',
        license='MIT',
        author='Asher Preska Steinberg',
        author_email='apsteinberg@nyu.edu',
        packages=['collectresults'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['getDivergence=collectresults.cli:main'],
            }
      )
