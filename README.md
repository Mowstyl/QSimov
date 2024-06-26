# QSimov

QSimov is a quantum computing toolkit developed in Python, initially as an End Degree Project at UCLM. It allows the user to work both at a low level using quantum gates and registers and at a higher level using the design tools. It has been created in such a way that it allows an almost immediate translation from the circuit diagram to the code necessary to execute it.
From the moment of its creation, QSimov has been focused on being able to be used in learning and not on big computers and datacenters. For this reason, the simplicity of use and the low hardware requirements have always been what have guided its development. It is capable of running up to 29 entangled qubits on a laptop with 16 GB of RAM and several thousand non-entangled qubits under the same conditions.

## Getting Started

By following this instructions you will have a working binary.

### Prerequisites

CPython>=3.8 x64 is needed to run QSimov with its default core, Doki.
You also need the following python libraries in order to use QSimov:
 - numpy>=1.21
 - matplotlib>=3.5.1
 - doki-Mowstyl>=1.5.2

### Building from sources

You need the source code, either by cloning the project or by downloading it. You also need the "wheels" and the "build" python modules.
Execute the following command while on the root directory of the project to generate the .whl file into the dist folder:
 - python -m build

### Installing

You may install QSimov with pip by downloading the latest version from PyPI with
 - pip install qsimov
You can also install the generated .whl package, in case you are building it from the source, with the following command.
 - pip install dist/qsimov-%QSimovVersion%-py3-none-any.whl

## Testing

The full test suite can be run with:
```
python tests/qtest.py <min> <max> [seed] [verbose]
```
- min: minimum number of qubits to use in the tests involving data structures with a variable number of qubits. min has to be at least 3.
- max: minimum number of qubits to use in the tests involving data structures with a variable number of qubits. Don't use more than 8 since the simulator tests against operations made in numpy, which is not made for quantum computing simulation and therefore can be very slow with a number of qubits greater than 7.
- seed (optional): the seed to use in random operations. If none specified, the used seed will be printed in console so you can always repeat a failed test.
- verbose (optional): whether to print extra information.

## Contributing

Feel free to create pull requests. Keep in mind that the requests must be synced to the *main* branch. If they have conflicts they'll be rejected.
Also, the pull requests must have a brief explaination of the work done.
They must be free of any binary file and the code has to be indented (plain code won't be accepted).

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags](https://github.com/Mowstyl/QSimovC/tags) section on this repository.

## Authors

* **Hernán Indíbil de la Cruz Calvo** - *Initial work* - [Mowstyl](https://github.com/Mowstyl)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

When a project has no LICENSE.md file, it's licensed under the most restrictive license (that means copyright). This is **always** true, so keep it in mind when using someone else's code.
At the moment, MIT License is being used.

## Acknowledgments

* **Billie Thompson** - *This README.md was based on the one provided by this person. Great work!* - [PurpleBooth](https://github.com/PurpleBooth)
* **StackEdit** - *Markdown editor used* - [StackEdit](https://stackedit.io/)
