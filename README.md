# Instructions:

1. Download strain_calc.exe from releases

2. Run using command line:

```
./strain_calc.exe <displacement file path> <subset size>
```

3. Program will output strain in a txt file

Note: displacement.txt file must be formatted in the following fashion:

```
<pixelX> <pixelY> <motionX> <motionY>

```
And must increment pixelY values before pixelX values

Example:

```
0 0 1.343 2.555
0 1 5.333 2.344
...
1 0 4.323 3.234
1 1 1.123 5.343 
...
```


## Building

### Dependencies

Eigen is used for matrix operations. It can be downloaded from the [Eigen website](http://eigen.tuxfamily.org/index.php).
Put the Eigen folder in the dependencies folder at `dependencies/Eigen/`.

### Building

To build the program, run the following command:

```bash
make
```

This will create a `strain_calc` executable in the current directory.



## Credits

- [Eigen](http://eigen.tuxfamily.org/index.php) for matrix operations
- [cnpy](https://github.com/rogersce/cnpy) for reading and writing numpy files