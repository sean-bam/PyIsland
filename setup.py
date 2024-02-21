import glob
from setuptools import find_packages, setup

requirements = [
    req.strip()
    for req in open("requirements.txt", "r").readlines()
    if not req.startswith("#")
]

setup(
    name="PyIsland",
    scripts=[script for script in glob.glob("bin/*")],
    packages=find_packages(),
    install_requires=requirements,
    version="0.0.1",
    description="Libary for comparative genomics",
    author="SMB",
    license="",
)
