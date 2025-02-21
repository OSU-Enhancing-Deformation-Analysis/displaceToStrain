// strain_calc.cpp : Defines the entry point for the application.
//

#include "strain_calc.h"
#include "dependencies/Eigen/Dense"
#include "iostream"
#include "matplot/matplot.h"

using namespace std;
using namespace Eigen;


/*
Process:
1. Split array of displacements into subsets
2. For each subset:
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

int getSubsets(vector<vector<pair<double, double>>> motionArray, int dims[2], int subsetSize, vector<pair<int, int>>& subsetLocs, vector<vector<pair<double, double>>>& displacements, vector<vector<pair<double, double>>>& distances);

int main()
{
    int dims[2] = { 100, 100 };
    std::vector<std::vector<std::pair<double, double>>> motionArray(
        dims[0],
        std::vector<std::pair<double, double>>(dims[1], pair(0.0, 0.0))
    );

    //------------fake motion for testing (thanks chatgpt): -------------
    // Seed random number generator for noise.
    srand(static_cast<unsigned int>(time(nullptr)));

    // Simulation parameters
    double strain = 0.005;      // 0.5% strain (tension) in the x direction
    double poisson = 0.3;       // Typical Poisson's ratio (contraction in y)
    double noiseLevel = 0.0001; // Small random noise level

    // Fill the motionArray with displacement values
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            // Assume that the pixel position (i, j) corresponds to (x, y)
            double x = static_cast<double>(i);
            double y = static_cast<double>(j);

            // Compute displacements:
            double u = strain * pow(x, 1.2) + pow(y,1.05);   // Now u depends on y
            double v = -poisson * strain * pow(y, 1.2) + pow(x,1.05); // And v depends on x

            // Optionally add some random noise to simulate measurement uncertainty:
            double noise_u = noiseLevel * ((rand() % 1000) / 1000.0 - 0.5);
            double noise_v = noiseLevel * ((rand() % 1000) / 1000.0 - 0.5);
            u += noise_u;
            v += noise_v;

            motionArray[i][j] = std::make_pair(u, v);
        }
    }
    //-----------------------------------------------------------------------


    int subsetSize = 5;
    vector<pair<int, int>> subsetLocs;
    vector<vector<pair<double, double>>> displacements;
    vector<vector<pair<double, double>>> distances;
    int numSubsets = getSubsets(motionArray, dims, subsetSize, subsetLocs, displacements, distances);
    int numNeigh = subsetSize * subsetSize - 1;

    /*for (int i = 0; i < subsetLocs.size(); i++) {
        cout << "x: " << subsetLocs[i].first << " y: " << subsetLocs[i].second << endl;
    }*/

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
        double coeffsX[3] = {0.0};
        double coeffsY[3] = {0.0};
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

    int subXSize = dims[0] / subsetSize;
	int subYSize = dims[1] / subsetSize;

    // Create motion arrays for plotting
    vector<vector<double>> motionX(dims[1], vector<double>(dims[0]));
    vector<vector<double>> motionY(dims[1], vector<double>(dims[0]));
    for (int x = 0; x < dims[0]; ++x) {
        for (int y = 0; y < dims[1]; ++y) {
            motionX[y][x] = motionArray[x][y].first;
            motionY[y][x] = motionArray[x][y].second;
        }
    }

	// Create strain Arrays
	vector<vector<double>> strainArrXX(subYSize, vector<double>(subXSize));
	vector<vector<double>> strainArrYY(subYSize, vector<double>(subXSize));
	vector<vector<double>> strainArrXY(subYSize, vector<double>(subXSize));

	for (int i = 0; i < numSubsets; ++i) {
		int x = subsetLocs[i].first;
		int y = subsetLocs[i].second;
		strainArrXX[y][x] = std::get<0>(strains[i]);
		strainArrYY[y][x] = std::get<1>(strains[i]);
		strainArrXY[y][x] = std::get<2>(strains[i]);
        //cout << "x " << std::get<0>(strains[i]) << "y " << std::get<1>(strains[i]) << "xy " << std::get<2>(strains[i]) << endl;
	}

    matplot::subplot(2, 3, 1);
    matplot::imagesc(motionX);
    matplot::title("Motion X");
    matplot::colorbar();

    matplot::subplot(2, 3, 2);
    matplot::imagesc(motionY);
    matplot::title("Motion Y");
    matplot::colorbar();

    matplot::subplot(2, 3, 3);
    matplot::imagesc(strainArrXX);
    matplot::title("Strain XX");
    matplot::colorbar();

    matplot::subplot(2, 3, 4);
    matplot::imagesc(strainArrYY);
    matplot::title("Strain YY");
    matplot::colorbar();

    matplot::subplot(2, 3, 5);
    matplot::imagesc(strainArrXY);
    matplot::title("Strain XY");
    matplot::colorbar();

    matplot::show();


	return 0;
}


//returns number of subsets
int getSubsets(vector<vector<pair<double, double>>> motionArray, int dims[2], int subsetSize, vector<pair<int, int>> &subsetLocs, vector<vector<pair<double, double>>> &displacements, vector<vector<pair<double, double>>> &distances) {
    int subsetsNum[2] = { dims[0] / subsetSize, dims[1] / subsetSize };
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