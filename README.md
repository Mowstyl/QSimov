# QSimov

This is the quantum framework QSimov, written in both Python and C.

## Getting Started

By following this instructions you will have a working binary.

### Prerequisites

You only need the source code, either by cloning the project or by downloading it, and Python >= 3.7. The QSimov quantum computing simulator works better when used under a x64 architecture. Expect lower maximum number of QuBits when on a x86 system.

### Installing

If you've tried before to build the project, you may want to clean it before trying to build it again.
If it's the first time you build it, it's safe to skip this step.
The cleaning can be performed by using
 - **Windows**
You don't need to do anything special, just extract the source code in a folder.
 - **GNU-Linux**
The lib folder has to be in the LDD search path. You can do that either by running
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/lib
```
before executing the simulator or, if that fails or if you want it to be always available, executing the *addLibraryPath.sh* script as admin.

## Testing

The full test suite can be run with:
```
python qtest.py <min> <max> [seed]
```
- min: minimum number of qubits to use in the tests involving data structures with a variable number of qubits. min has to be at least 3.
- max: minimum number of qubits to use in the tests involving data structures with a variable number of qubits. Don't use more than 8 since the simulator tests against operations made in numpy, which is not made for quantum computing simulation and therefore can be very slow with a number of qubits greater than 7.
- seed (optional): the seed to use in random operations. If none specified, the used seed will be printed in console so you can always repeat a failed test.

## libqsimov and libfunmant built With

* [MSYS2](https://www.msys2.org/) - Software distro and building platform for Windows
* MinGW-w64 - Build tools for Windows systems. Obtained with MSYS2 to get the latest gcc version.
* [GCC](https://gcc.gnu.org/) - C compiler for GNU-Linux systems.

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
