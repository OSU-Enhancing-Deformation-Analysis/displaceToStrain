// strain_calc.cpp : Defines the entry point for the application.
//
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "strain_calc.h"
#include "Eigen/Dense"
#include "cnpy.h"
#include "stb_image_write.h"
#include "iostream"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <regex>
#include <vector>
#include <string>
#include <tuple>

const int TXT = 0;
const int NPY = 1;
const int JPG = 2;

using namespace std;
using namespace Eigen;
namespace fs = filesystem;

void processFile(const string& filePath, const string& outputDir, int subsetSize, int exportFormat);
string getFileName(const string& filePath);
vector<vector<pair<double, double>>> getMotionFromFile(string fileName);
vector<vector<pair<double, double>>> getMotionNumpyFromFile(string fileName);
vector<vector<pair<double, double>>> getMotionFromTiles(string folderName, int tileXNum, int tileYNum, int tileW, int tileH);
pair<int, int> extractTileCoordinates(const string& filename);
void exportStrainToFile(vector<vector<tuple<double, double, double>>> strainArray, string fileName);
void exportStrainToNumpy(vector<vector<tuple<double, double, double>>> strainArray, string fileName);
void exportStraintoImage(vector<vector<tuple<double, double, double>>> strainArray, string fileName);
void exportStraintoImages(vector<vector<tuple<double, double, double>>> strainArray, string fileNameXX, string fileNameYY, string fileNameXY);
pair<double, double> minMax(vector<vector<tuple<double, double, double>>> strainArray);
tuple<double, double, double> absMax(vector<vector<tuple<double, double, double>>> strainArray);
static double normalize(double min, double max, double val);
int getSubsets(vector<vector<pair<double, double>>> motionArray, int subsetSize, vector<pair<int, int>>& subsetLocs, vector<vector<pair<double, double>>>& displacements, vector<vector<pair<double, double>>>& distances);
vector<vector<tuple<double, double, double>>> calcStrains(vector<vector<pair<double, double>>> motionArray, int subsetSize);

//arguments: <disp file name> <subset size> [-npy] [-o <output_dir>]
int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Not enough arguments." << endl;
        cout << "Format:" << endl;
        cout << "./strain_calc.exe <disp file/folder path> <subset size> [-npy/jpg] [-o <output_dir>]" << endl;
        return 0;
    }
    string filePath = argv[1];
    int subsetSize = stoi(argv[2]);
    int exportFormat = TXT;
    string outputDir = ""; // Default to current directory

    int arg_index = 3;
    while (arg_index < argc) {
        string arg = argv[arg_index];
        if (arg == "-npy" || arg == "--numpy") {
            exportFormat = NPY;
        } else if (arg == "-jpg") {
            exportFormat = JPG;
        } else if (arg == "-o" || arg == "--output_dir") {
            if (arg_index + 1 < argc) {
                outputDir = argv[arg_index + 1];
                arg_index += 1; // Skip the next argument as it's the directory
                // Ensure outputDir ends with a directory separator
                if (!outputDir.empty() && outputDir.back() != '/' && outputDir.back() != '\\') {
                    outputDir += "/";
                }
            } else {
                cerr << "Error: -o or --output_dir option requires a directory path." << endl;
                return 1;
            }
        }
        arg_index += 1;
    }

    if (fs::is_regular_file(filePath)) { //single file
        cout << "Processing single file" << endl;
        processFile(filePath, outputDir, subsetSize, exportFormat);
    }
    else if (fs::is_directory(filePath)) { //folder
        //create output folder
        fs::path pathObj(filePath);
        outputDir = outputDir + pathObj.filename().string() + "_output/";
        fs::create_directory(outputDir);
        //go through each file in input folder
        cout << "Running on all files in folder... (may take a while)" << endl;
        for (const auto& entry : fs::directory_iterator(filePath)) {
            if (entry.is_regular_file()) {
                string fileName = entry.path().filename().string();
                string extension = fileName.substr(fileName.find_last_of(".") + 1);
                if (extension != "txt" && extension != "npy") {
                    continue;
                }
                processFile(entry.path().string(), outputDir, subsetSize, exportFormat);
            }
        }
    }


}

void processFile(const string& filePath, const string& outputDir, int subsetSize, int exportFormat)
{
    // Check if the file ends with .npy
    vector<vector<pair<double, double>>> motionArray;
    if (filePath.substr(filePath.size() - 4) == ".npy") {

        motionArray = getMotionNumpyFromFile(filePath);
    }
    else {
        motionArray = getMotionFromFile(filePath);
    }

    vector<vector<tuple<double, double, double>>> strainArray = calcStrains(motionArray, subsetSize);

    string baseFileName = getFileName(filePath) + "_strain_result";

    switch (exportFormat) {
    case TXT:
        exportStrainToFile(strainArray, outputDir + baseFileName + ".txt");
        break;
    case NPY:
        exportStrainToNumpy(strainArray, outputDir + baseFileName + ".npy");
        break;
    case JPG:
        string filePathXX = outputDir + "XX_" + baseFileName + ".jpg";
        string filePathYY = outputDir + "YY_" + baseFileName + ".jpg";
        string filePathXY = outputDir + "XY_" + baseFileName + ".jpg";
        exportStraintoImages(strainArray, filePathXX, filePathYY, filePathXY);
        break;
    }
}

string getFileName(const string& filePath) {
    // Find the last occurrence of directory separator
    size_t lastSlash = filePath.find_last_of("/\\");

    // Extract the file name from the path
    string fileName = (lastSlash == string::npos) ? filePath : filePath.substr(lastSlash + 1);

    // Find the last occurrence of a dot (.)
    size_t lastDot = fileName.find_last_of(".");

    // Extract only the file name without the extension
    if (lastDot != string::npos) {
        fileName = fileName.substr(0, lastDot);
    }

    return fileName;
}

vector<vector<pair<double, double>>> getMotionFromFile(string fileName)
{
    // Open the  file
    ifstream file(fileName);

    // Check if file is open
    if (!file.is_open()) {
        cerr << "Error opening file." << endl;
        return {};
    }

    string line;

    // Read the file line by line
    vector<vector<pair<double, double>>> motionArray;
    while (getline(file, line)) {
        stringstream ss(line);  // Create a stringstream from the line
        string value;
        vector<string> row;

        // Read each value (separated by commas) and store in a vector
        while (getline(ss, value, ' ')) {
            row.push_back(value);
        }

        //add row values to motionArray based on x and y value
        if (row.size() != 4) { //invalid row
            continue;
        }
        int x = stoi(row[0]);
        int y = stoi(row[1]);
        float motionX = stof(row[2]);
        float motionY = stof(row[3]);

        if (x >= motionArray.size()) {
            motionArray.push_back(vector<pair<double, double>>());
        }

        motionArray[x].push_back(pair(motionX, motionY));

    }

    // Close the file
    file.close();

    return motionArray;
}

vector<vector<pair<double, double>>> getMotionNumpyFromFile(string fileName) {
	cnpy::NpyArray arr = cnpy::npy_load(fileName);

	if (arr.shape.size() != 3 || arr.shape[2] != 2) {
		cerr << "Invalid .npy file format. Expected shape (rows, cols, 2)." << endl;
		return {};
	}

	// Check if the data type matches double
	bool isDouble = true;
	if (arr.word_size != sizeof(double)) {
		isDouble = false;
	}

	size_t rows = arr.shape[0];
	size_t cols = arr.shape[1];

	vector<vector<pair<double, double>>> motionArray(rows, vector<pair<double, double>>(cols));

	// Read the data assuming row-major order.
    if (isDouble) {
        double* data = arr.data<double>();
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                size_t index = (i * cols + j) * 2;  // Each element is 2 doubles (x, y)
                motionArray[i][j] = { data[index], data[index + 1] };
            }
        }
    }
    else {
        float* data = arr.data<float>();
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                size_t index = (i * cols + j) * 2;  // Each element is 2 floats (x, y)
                motionArray[i][j] = { double(data[index]), double(data[index + 1]) };
            }
        }
    }
	
	return motionArray;
}

vector<vector<pair<double, double>>> getMotionFromTiles(string folderName, int tileXNum, int tileYNum, int tileW, int tileH)
{
    vector<vector<pair<double, double>>> motionArray(tileXNum, vector<pair<double, double>>(tileYNum, pair(0.0, 0.0)));

    // Iterate over all files in the folder
    for (const auto& entry : fs::directory_iterator(folderName)) {
        if (entry.is_regular_file()) {
            ifstream file(entry.path());

            // Check if file is open
            if (!file.is_open()) {
                cerr << "Error opening file: " << entry.path() << endl;
                continue;
            }

            string line;

            // Read the file line by line
            while (getline(file, line)) {
                stringstream ss(line);  // Create a stringstream from the line
                string value;
                vector<string> row;

                // Read each value (separated by space) and store in a vector
                while (getline(ss, value, ' ')) {
                    row.push_back(value);
                }

                // Skip invalid rows
                if (row.size() != 4) {
                    continue;
                }

                int x = stoi(row[0]);
                int y = stoi(row[1]);
                float motionX = stof(row[2]);
                float motionY = stof(row[3]);

                // Adjust x and y coordinates based on the tile position
                pair<int, int> tileCoords = extractTileCoordinates(entry.path().filename().string());
                int tileIndexX = tileCoords.first;
                int tileIndexY = tileCoords.second;

                // Use tileIndexX and tileIndexY for adjusting global coordinates
                int globalX = x + tileIndexX * tileW;
                int globalY = y + tileIndexY * tileH;

                motionArray[globalX][globalY] = pair(motionX, motionY);

                // Add the motion data to the correct position in the larger array
                motionArray[globalX][globalY] = pair(motionX, motionY);
            }

            // Close the file
            file.close();
        }
    }

    return motionArray;
}

pair<int, int> extractTileCoordinates(const string& filename) {
    // Define the regex pattern to match the suffix _x_y
    regex pattern("_(\\d+)_(\\d+)\\.txt$");

    smatch matches;
    if (regex_search(filename, matches, pattern)) {
        // Extract the x and y coordinates from the match
        int x = stoi(matches[1].str());
        int y = stoi(matches[2].str());
        return { x, y };
    }

    // If the filename doesn't match the expected pattern, return -1, -1 as an error code
    return { -1, -1 };
}

void exportStrainToFile(vector<vector<tuple<double, double, double>>> strainArray, string fileName)
{
    ofstream file(fileName);

    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    for (int y = 0; y < strainArray[0].size(); ++y) {
        for (int x = 0; x < strainArray.size(); ++x) {
            file << x << " " << y << " " <<
                get<0>(strainArray[x][y]) << " " <<
                get<1>(strainArray[x][y]) << " " <<
                get<2>(strainArray[x][y]) << endl;
        }
    }

    file.close();
    cout << "txt file exported successfully!\n";
}

void exportStrainToNumpy(vector<vector<tuple<double, double, double>>> strainArray, string fileName) {
    size_t rows = strainArray.size();
    size_t cols = (rows > 0) ? strainArray[0].size() : 0; // Handle empty array case

    if (rows == 0 || cols == 0) {
        cerr << "Warning: Strain array is empty, nothing to export to numpy." << endl;
        return;
    }

    vector<double> numpyData;
    for (size_t x = 0; x < rows; ++x) {
        for (size_t y = 0; y < cols; ++y) {
            numpyData.push_back(get<0>(strainArray[x][y]));
            numpyData.push_back(get<1>(strainArray[x][y]));
            numpyData.push_back(get<2>(strainArray[x][y]));
        }
    }

    vector<size_t> shape = {rows, cols, 3};
	cnpy::npy_save(fileName, &numpyData[0], shape, "w");
	cout << ".npy file exported successfully with shape " << shape[0] << "x" << shape[1] << "x" << shape[2] << endl;
}

void exportStraintoImage(vector<vector<tuple<double, double, double>>> strainArray, string fileName)
{
	int width = strainArray[0].size(), height = strainArray.size();
	unsigned char* image_data = new unsigned char[width * height * 3];  // RGB image

    pair<double, double> minMaxVals = minMax(strainArray);

	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
            int index = (width * y + x) * 3;
			image_data[index + 0] = static_cast<unsigned char>(normalize(minMaxVals.first, minMaxVals.second, get<0>(strainArray[y][x])) * 255);
			image_data[index + 1] = static_cast<unsigned char>(normalize(minMaxVals.first, minMaxVals.second, get<1>(strainArray[y][x])) * 255);
			image_data[index + 2] = static_cast<unsigned char>(normalize(minMaxVals.first, minMaxVals.second, get<2>(strainArray[y][x])) * 255);
		}
	}


    // Fill image_data with pixel values
    stbi_write_jpg(fileName.c_str(), width, height, 3, image_data, 100);
    delete[] image_data;
}

void exportStraintoImages(vector<vector<tuple<double, double, double>>> strainArray, string fileNameXX, string fileNameYY, string fileNameXY)
{
    int width = strainArray[0].size(), height = strainArray.size();
    unsigned char* imageDataXX = new unsigned char[width * height * 3];  // RGB image
    unsigned char* imageDataYY = new unsigned char[width * height * 3];  // RGB image
    unsigned char* imageDataXY = new unsigned char[width * height * 3];  // RGB image

    tuple<double, double, double> maxVals = absMax(strainArray);

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            int index = (width * y + x) * 3;
            imageDataXX[index + 0] = 0;
            imageDataXX[index + 1] = 0;
            imageDataXX[index + 2] = 0;
            imageDataYY[index + 0] = 0;
            imageDataYY[index + 1] = 0;
            imageDataYY[index + 2] = 0;
            imageDataXY[index + 0] = 0;
            imageDataXY[index + 1] = 0;
            imageDataXY[index + 2] = 0;
            if (get<0>(strainArray[y][x]) < 0) {
                imageDataXX[index + 0] = abs(get<0>(strainArray[y][x])) / get<0>(maxVals) * 255;
            }
            else {
                imageDataXX[index + 2] = get<0>(strainArray[y][x]) / get<0>(maxVals) * 255;
            }
            if (get<1>(strainArray[y][x]) < 0) {
                imageDataYY[index + 0] = abs(get<1>(strainArray[y][x])) / get<1>(maxVals) * 255;
            }
            else {
                imageDataYY[index + 2] = get<1>(strainArray[y][x]) / get<1>(maxVals) * 255;
            }
            if (get<2>(strainArray[y][x]) < 0) {
                imageDataXY[index + 0] = abs(get<2>(strainArray[y][x])) / get<2>(maxVals) * 255;
            }
            else {
                imageDataXY[index + 2] = get<2>(strainArray[y][x]) / get<2>(maxVals) * 255;
            }
        }
    }


    // Fill image_data with pixel values
    stbi_write_jpg(fileNameXX.c_str(), width, height, 3, imageDataYY, 100);
    stbi_write_jpg(fileNameYY.c_str(), width, height, 3, imageDataXX, 100);
    stbi_write_jpg(fileNameXY.c_str(), width, height, 3, imageDataXY, 100);
    delete[] imageDataXX;
    delete[] imageDataYY;
    delete[] imageDataXY;
}

pair<double, double> minMax(vector<vector<tuple<double, double, double>>> strainArray) {
    double minVal = 10000;
    double maxVal = -10000;
    int width = strainArray.size(), height = strainArray[1].size();
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            minVal = get<0>(strainArray[x][y]) < minVal ? get<0>(strainArray[x][y]) : minVal;
            minVal = get<1>(strainArray[x][y]) < minVal ? get<1>(strainArray[x][y]) : minVal;
            minVal = get<2>(strainArray[x][y]) < minVal ? get<2>(strainArray[x][y]) : minVal;
            maxVal = get<0>(strainArray[x][y]) > maxVal ? get<0>(strainArray[x][y]) : maxVal;
            maxVal = get<1>(strainArray[x][y]) > maxVal ? get<1>(strainArray[x][y]) : maxVal;
            maxVal = get<2>(strainArray[x][y]) > maxVal ? get<2>(strainArray[x][y]) : maxVal;
        }
    }
    return pair<double, double>(minVal, maxVal);
}
tuple<double, double, double> absMax(vector<vector<tuple<double, double, double>>> strainArray) {
    tuple<double, double, double> maxVals = tuple<double, double, double>(0, 0, 0);
    int width = strainArray.size(), height = strainArray[1].size();
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            get<0>(maxVals) = get<0>(maxVals) < abs(get<0>(strainArray[x][y])) ? abs(get<0>(strainArray[x][y])) : get<0>(maxVals);
            get<1>(maxVals) = get<1>(maxVals) < abs(get<1>(strainArray[x][y])) ? abs(get<1>(strainArray[x][y])) : get<1>(maxVals);
            get<2>(maxVals) = get<2>(maxVals) < abs(get<2>(strainArray[x][y])) ? abs(get<2>(strainArray[x][y])) : get<2>(maxVals);
        }
    }
    return maxVals;
}


static double normalize(double min, double max, double val) {
    return (val - min) / (max - min);
}


/*
motionArray: 2D array of motion vector pairs for each pixel
subsetSize: size of each (square) subset
subsetLocs: vector to store subset locations
displacements: vector to store list of displacements for each subset
distances: vector to store list of distances for each subset

returns number of subsets

This function seperates per-pixel motion data in distinct subsets for strain calculation.
May be updated to check for invalid or empty pixels later.
*/
int getSubsets(vector<vector<pair<double, double>>> motionArray, int subsetSize, vector<pair<int, int>> &subsetLocs, vector<vector<pair<double, double>>> &displacements, vector<vector<pair<double, double>>> &distances) {
    int sizeX = motionArray.size();
    int sizeY = motionArray[0].size();
    //num subsets rounded down if image is not divisible by subset
    //this leaves a small portion of the image unprocessed
    int subsetsNum[2] = { int(sizeX / subsetSize), int(sizeY / subsetSize) };
    for (int x = 0; x < subsetsNum[0]; x++) {
        for (int y = 0; y < subsetsNum[1]; y++) {
            int xCenter = x * subsetSize + int(subsetSize / 2);
            int yCenter = y * subsetSize + int(subsetSize / 2);
            subsetLocs.push_back(pair(x, y));
            vector<pair<double, double>> subDisplacements;
            vector<pair<double, double>> subDistances;
            for (int px = x * subsetSize; px < (x + 1) * subsetSize; px++) {
                for (int py = y * subsetSize; py < (y + 1) * subsetSize; py++) {
                    if (px == xCenter and py == yCenter) {
                        continue;
                    }
                    subDisplacements.push_back(pair(motionArray[px][py].first, motionArray[px][py].second));
                    int xDist = px - xCenter;
                    int yDist = py - yCenter;
                    subDistances.push_back(pair(xDist, yDist));
                }
            }
            displacements.push_back(subDisplacements);
            distances.push_back(subDistances);
        }
    }
    return subsetsNum[0] * subsetsNum[1];
}

/*
motionArray: 2D array of motion vector pairs for each pixel
subsetSize: size of each (square) subset

returns 2D array of tuples containing (strainXX, strainYY, strainXY)
(this array is of dimension (numSubsetsX, numSubsetsY)

This function calculates strain values for each subset.
The general process is as follows:
For each subset:
    1. Solve least squares fit for displacement given neighbor distance
        X^T = [1.0, dist_x, dist_y]
        u_x = [disp_x]
        u_y = [disp_y]
        coeffs_x = [dudx, dudy] = (X^T * X)^-1 * X^T * u_x
        coeffs_y = [dvdx, dvdy] = (X^T * X)^-1 * X^T * u_y
    2. compute green lagrange strain
        strain_xx = 0.5*(2.0*dudx + dudx*dudx + dvdx*dvdx);
        strain_yy = 0.5*(2.0*dvdy + dudy*dudy + dvdy*dvdy);
        strain_xy = 0.5*(dudy + dvdx + dudx*dudy + dvdx*dvdy);
*/
vector<vector<tuple<double, double, double>>> calcStrains(vector<vector<pair<double, double>>> motionArray, int subsetSize)
{
    vector<pair<int, int>> subsetLocs;
    vector<vector<pair<double, double>>> displacements;
    vector<vector<pair<double, double>>> distances;
    int numSubsets = getSubsets(motionArray, subsetSize, subsetLocs, displacements, distances);
    int numNeigh = subsetSize * subsetSize - 1;

    vector<tuple<double, double, double>> strains;

    for (int subset = 0; subset < numSubsets; subset++) {
        //for each subset
        //calculate X^T matrix and ux, uy
        MatrixXd  Xt = MatrixXd::Zero(3, numNeigh);
        VectorXd  ux = VectorXd::Zero(numNeigh);
        VectorXd  uy = VectorXd::Zero(numNeigh);
        for (int neigh = 0; neigh < numNeigh; neigh++) {
            //displacements
            ux(neigh) = displacements[subset][neigh].first;
            uy(neigh) = displacements[subset][neigh].second;
            //distance matrix
            Xt(0, neigh) = 1.0;
            Xt(1, neigh) = distances[subset][neigh].first;
            Xt(2, neigh) = distances[subset][neigh].second;
        }
        //calculate X^T*X
        MatrixXd XtX = MatrixXd::Zero(3, 3);
        for (int k = 0; k < 3; ++k) {
            for (int m = 0; m < 3; ++m) {
                for (int neigh = 0; neigh < numNeigh; ++neigh) {
                    XtX(k, m) += Xt(k, neigh) * Xt(m, neigh);
                }
            }
        }

        //calculate inverse
        MatrixXd XtXi = XtX.inverse();

        //compute X^T*u
        VectorXd Xtux = VectorXd::Zero(3);
        VectorXd Xtuy = VectorXd::Zero(3);
        for (int i = 0; i < 3; ++i) {
            for (int neigh = 0; neigh < numNeigh; ++neigh) {
                Xtux(i) += Xt(i, neigh) * ux(neigh);
                Xtuy(i) += Xt(i, neigh) * uy(neigh);
            }
        }

        //calculate coefficients
        double coeffsX[3] = { 0.0 };
        double coeffsY[3] = { 0.0 };
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                coeffsX[i] += XtXi(i, j) * Xtux(j);
                coeffsY[i] += XtXi(i, j) * Xtuy(j);
            }
        }

        //calculate strain
        double dudx = coeffsX[1];
        double dudy = coeffsX[2];
        double dvdx = coeffsY[1];
        double dvdy = coeffsY[2];

        double strainXX = 0.5 * (2.0 * dudx + dudx * dudx + dvdx * dvdx);
        double strainYY = 0.5 * (2.0 * dvdy + dudy * dudy + dvdy * dvdy);
        double strainXY = 0.5 * (dudy + dvdx + dudx * dudy + dvdx * dvdy);

        strains.push_back(tuple(strainXX, strainYY, strainXY));
    }

    int sizeX = motionArray.size();
    int sizeY = motionArray[0].size();
    int numSubsetsX = sizeX / subsetSize;
    int numSubsetsY = sizeY / subsetSize;

    // Create 2D strain Array
    vector<vector<tuple<double, double, double>>> strainArr(numSubsetsX, vector<tuple<double, double, double>>(numSubsetsY, tuple(0.0, 0.0, 0.0))); 

    for (int i = 0; i < numSubsets; ++i) {
        int x = subsetLocs[i].first;
        int y = subsetLocs[i].second;
        get<0>(strainArr[x][y]) = get<0>(strains[i]);
        get<1>(strainArr[x][y]) = get<1>(strains[i]);
        get<2>(strainArr[x][y]) = get<2>(strains[i]);
    }
    return strainArr;
}