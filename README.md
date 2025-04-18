# Instructions:

1. Download strain_calc.exe from releases

2. Run using command line:

```
./strain_calc.exe <disp file/folder path> <subset size> [-npy/jpg] [-o <output_dir>]
```

3. Program will output strain in a txt file by default, or you can specify -npy or -jpg to output in other formats.

> Note: displacement.txt file must be formatted in the following fashion:

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

> Note: When using the -npy flag, the displacement file must be a .npy file. The values should be double precision and in the shape (rows, cols, 2).



## Building

### Dependencies

Eigen is used for matrix operations. It can be downloaded from the [Eigen website](http://eigen.tuxfamily.org/index.php).
Put the Eigen folder in the dependencies folder at `dependencies/Eigen/`.

OpenCV is required for strain pre/post processing. It can be downloaded from [OpenCV website](https://opencv.org/releases/).
**After downloading, you must:**

1. Add the path to opencv/build on your computer to your PATH environment variable

2. In opencv/build/x64/vc16/bin, copy all dll and pdb files to the same folder as your strain_calc.exe

### Building

To build the program, run the following command:

```bash
make
```

This will create a `strain_calc` executable in the current directory.



## Credits

- [Eigen](http://eigen.tuxfamily.org/index.php) for matrix operations
- [cnpy](https://github.com/rogersce/cnpy) for reading and writing numpy files
- [stb](https://github.com/nothings/stb/tree/master) for exporting images
- [OpenCV](https://github.com/opencv/opencv/tree/4.11.0) for pre/post processing
