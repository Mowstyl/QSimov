"""Installation module."""
from setuptools import setup, find_packages


def main():
    """Code to be executed on install."""
    setup(
        name="qsimov-Mowstyl",
        version="4.3.0",
        author="Hernán Indíbil de la Cruz Calvo",
        author_email="HernanIndibil.LaCruz@alu.uclm.es",
        license="MIT",
        packages=find_packages(include=['qsimov', 'qsimov.*']),
        url="https://github.com/Mowstyl/QSimov",
        description="QSimov Quantum computer simulator",
        long_description="QSimov is a quantum computer simulator based on " +
                         "the circuit model.",
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Software Development :: Libraries :: Python Modules',
        ],
        keywords="quantum",
        install_requires=[
            "numpy>=1.21",
            "matplotlib>=3.5.1",
            "doki-Mowstyl>=1.3.2"
        ]
    )


if __name__ == "__main__":
    main()
