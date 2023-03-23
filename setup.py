from setuptools import find_packages, setup
setup(
    name='Paragon',
    packages=find_packages(include=['Paragon']),
    version='0.1.0',
    description='Network Inference Lib',
    author='M Kaan ARICI, Nurcan Tuncbag',
    license='MIT',
    url='https://github.com/metunetlab/Paragon',
    python_requires='>=3.6',
    install_requires=[
            "numpy==1.19.5",
            "networkx==2.5",
            "pandas==1.1.5",
            "community== 1.0.0b1",
            "scipy==1.5.4",
            ]
)
