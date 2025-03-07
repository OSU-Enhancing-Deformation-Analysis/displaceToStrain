Instructions:

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
