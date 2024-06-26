## setup.py file for vqt
from pathlib import Path

from distutils.core import setup

## Get the parent directory of this file
dir_parent = Path(__file__).parent

## Get requirements from requirements.txt
def read_requirements():
    with open(str(dir_parent / "requirements.txt"), "r") as req:
        content = req.read()  ## read the file
        requirements = content.split("\n") ## make a list of requirements split by (\n) which is the new line character

    ## Filter out any empty strings from the list
    requirements = [req for req in requirements if req]
    ## Filter out any lines starting with #
    requirements = [req for req in requirements if not req.startswith("#")]
    ## Remove any commas, quotation marks, and spaces from each requirement
    requirements = [req.replace(",", "").replace("\"", "").replace("\'", "").strip() for req in requirements]

    return requirements
deps_all = read_requirements()


setup(
    name='Q_learning_2ABT',
    version=0.1,
    author='Richard Hakim, Daniel Hochbaum',
    license='LICENSE',

    packages=[
        'Q_learning_2ABT',
    ],
    
    install_requires=deps_all,
    # extras_require={},
)