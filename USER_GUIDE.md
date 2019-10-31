# QSimov

This is the quantum framework QSimov, written in both Python and C.

## Getting Started

By following this instructions you will have a working binary.

### Prerequisites

You only need the source code, either by cloning the project or by downloading it, and a C compiler. The QSimov quantum computing simulator works better when used under a x64 architecture. Expect lower maximum number of QuBits when on a x86 system.

 - **Windows**
To compile the project in Windows (without any changes to the makefile) you'll need a compiler that supports C99. The project has been tested to work with MinGW-64. The lastest version of Cygwin (last checked 20/11/2018) has an older version of GCC that doesn't support C18, but it still works with C99.

 - **GNU-Linux**
To compile the project you'll need a C compiler that supports C99 (again). GCC versions 4 has been tested to work with the provided makefile (last checked 20/11/2018).

### Customizing

You can change the data types used on the calculations performed by QSimov, as well as the print format and some other little customizations. In order to do that, you can edit the *headers/precision.h* file. Read the comments and change the code having them in mind. Choose wisely. Using *long double* precision has some unexpected results in some systems, so be sure to run some tests to check they work in your system.

### Installing

If you've tried before to build the project, you may want to clean it before trying to build it again.
If it's the first time you build it, it's safe to skip this step.
The cleaning can be performed by using
 - **Windows** (with MinGW)
```
mingw32-make clean
```
 - **GNU-Linux** (with GCC)
```
make clean
```
After cleaning the project, you can safely build it.
 - **Windows** (with MinGW)
```
mingw32-make
```
 - **GNU-Linux** (with GCC)
```
make
```
You will have the binary in *bin/QSimovC.exe*.

## Testing

TBA

## Deployment

TBA

## Built With

* [MinGW-w64](https://sourceforge.net/projects/mingw-w64/) - Build tools for Windows systems
* [GCC](https://gcc.gnu.org/) - C compiler for GNU-Linux systems.

## Contributing

Feel free to create pull requests. Keep in mind that the requests must be synced to the *main* branch. If they have conflicts they'll be rejected.
Also, the pull requests must have a brief explaination of the work done.
They must be free of any binary file and the code has to be indented (plain code won't be accepted).
The code must also follow the GNU C coding style [guideline](https://developer.gnome.org/programming-guidelines/stable/c-coding-style.html.en).

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags](https://github.com/Mowstyl/QSimovC/tags) section on this repository.

## Authors

* **Hernán Indíbil de la Cruz Calvo** - *Initial work* - [Mowstyl](https://github.com/Mowstyl)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

When a project has no LICENSE.md file, it's licensed under the most restrictive license (that means copyright). This is **always** true, so keep it in mind when using someone else's code.
At the moment, Copyright is being used. In the future the license may change, maybe MIT license, maybe GNU, who knows?

## Acknowledgments

* **Billie Thompson** - *This README.md was based on the one provided by this person. Great work!* - [PurpleBooth](https://github.com/PurpleBooth)
* **StackEdit** - *Markdown editor used* - [StackEdit](https://stackedit.io/)
