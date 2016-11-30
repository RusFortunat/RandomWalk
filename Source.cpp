#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <random>
#include <limits>

using namespace std;

void shiftright(int myarray[], int myarray2[], int element, int size);
void shiftleft(int myarray[], int myarray2[], int element, int size);
void PressKeyToContinue();

// The whole purpose of the program is to reproduce the central limit theorem and plot the gaussian
int main() {
	string outFile("D:\\Science\\Research\\Random_Walk\\RandWalk_1d.txt"); // File for storing the output.
	ofstream out_stream;
	out_stream.open(outFile.c_str());

	if (out_stream.fail()) // Just in case
	{
		cout << "Error occurred during opening file!" << endl;
		return 1;
	}

	if (out_stream)
	{
		//// Greeting ////
		cout << "Hi there!\nLets simulate random walk in 1d lattice!" << endl << endl;
		//// Declaration of the variables ////
		int L, Lat, R, dt, N, nwalks, iwalk, istep, iter, bc;
		int a, b, c, d, e, n, m, Rest;
		double p, plain, x, xsum, xsqsum, r, rho;
		//// Specifying the parameters of the simulation ////
		cout << "First, enter parameters for your 1d lattice:" << endl;
		cout << "Size of the lattice: "; // The size of one dimensional array
		cin >> L;
		cout << "\nSet the probability for atom to go the right: ";
		cin >> p; // Sets the probability to go to the left (and to the rigth: 1 - p )
		cout << "\nHow many steps for each walk: ";
		cin >> dt;
		if (dt > L / 2) {
			cout << "Error! The number of time steps cannot be greater that the half of the lattice!!!";
			PressKeyToContinue();
			return 1;
		}
		cout << "\nHow many times you want to simulate the walk: ";
		cin >> nwalks;

		//// Random walk algorithm ////
		// Mersenne Twister RNG //
		random_device rd{};
		mt19937 RNG{ rd() };
		uniform_int_distribution<int> Latt(0, L - 1);
		uniform_real_distribution<double> Rand{ 0.0, 1.0 };
		xsum = 0; xsqsum = 0; N = 1;
		int* Lattice = new int[L]; // Our one dimensional lattice
		int* Counter = new int[L]; // Here we will store info how many times each cell were occupied during the single run
		double* Prob = new double[L]; // Array with probabilities
		for (n = 0; n < L; n++) // Setting counter to zero
		{
			Prob[n] = 0; // 
		}

		//// Beginning of the simulation ////
		for (iwalk = 0; iwalk < nwalks; iwalk++) { // Number of runs
			for (n = 0; n < L; n++)// Setting all lattice cells to be empty first
			{
				Lattice[n] = 0;
			}
			for (n = 0; n < L; n++) // Setting counter to zero
			{
				Counter[n] = 0; // 
			}
			Lattice[L / 2] = 1;
			Counter[L / 2] = 1;
			cout << "\nInitial positions of the attoms in the lattice:\n";
			for (n = 0; n < L; n++) {
				cout << Lattice[n];
			}cout << endl;

			// Beginning of the single run
			for (istep = 0; istep < dt; istep++) { // Number of steps for each run
					// Picking the random particle in the array
				for (n = 0; n < L; n++) {
					if (Lattice[n] != 0) { // Working only with non-empty cells
					// Throwing the due. Generates a random number between 0 and 1. If number is less than probability p, go to the right, otherwise go to the left
						r = Rand(RNG);
						if (r > p) {
							shiftleft(Lattice, Counter, n, L);
						}
						else {
							shiftright(Lattice, Counter, n, L);
						}
						n += L;
					}
				}
				// Displaying change of the positions
				for (n = 0; n < L; n++) {
					cout << Lattice[n];
				}cout << endl;
			}
			// End of the run

			//Computing output parameters
			// Probability distribution
			for (n = 0; n < L; n++) {
				Prob[n] += Counter[n];
			}

			// Displacement
			for (n = 0; n < L; n++) {
				if (Lattice[n] != 0) {
					x = n - L / 2;
					xsum += x; // Accumulates all final displacements to compute <x>
					xsqsum += x*x; // Accumulates squared displacement to compute MSD
				}
			}

			// Check if we lost/acquire particle(s) during the run
			Rest = 0;
			for (n = 0; n < L; n++) {
				Rest += Lattice[n];
			}
			cout << "\nTotal particles number: before " << N << ",  after  " << Rest << endl;
			if (N != Rest) {
				cout << "Partice lost/acquired during the run!!!" << endl;
				PressKeyToContinue();
				return 1;
			}
		}
		// End of the simulation

		//// Computing general output ////
		// Probability distribution
		for (n = 0; n < L; n++) {
			Prob[n] = Prob[n] / (nwalks*dt*N);
		}

		plain = 1.0 / (L); // Plain probability distribution
		if (51 < L&&L < 501) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
				n += 5;
			}
		}
		if (501 < L&&L < 5001) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
				n += 50;
			}
		}
		if (5001 < L&&L < 50000) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
				n += 500;
			}
		}
		cout << "\nCompare with the flat distribution 1/L: " << plain << endl << endl;
		cout << "\nAverage displacement: " << xsum / nwalks << endl;
		cout << endl << "\nMean squared displacement: " << xsqsum / nwalks << endl;
		cout << "\nCompare with the number of timesteps: " << dt << endl;

		//// Exporting the output to the text file ////
		// Header
		out_stream << "This is a document that stores the output from Random Walk program. The purpose of this program is to reproduce the central limit theorem.\n\nThe simulation input parameters:\n";
		// Simulation parameters
		out_stream << "Lattice size: " << L << " cells" << "\nNumber of atoms: " << N << "\nProbability to go left: " << p << "\nNumber of walks: " << nwalks << "\nNumber of steps in each walk: " << dt << "\nMean displacement: " << xsum / nwalks << "\nMean square displacement" << xsqsum / nwalks << ".\n\n Cell number | Probability to find an atom\n";
		// Output
		for (n = 0; n < L; n++) {
			out_stream << n - L / 2 << "		" << Prob[n] << endl;
		}
		out_stream.close();
	}

	PressKeyToContinue();
	return 0;
}


void shiftleft(int myarray[], int myarray2[], int element, int size) // Shifts atom to the left
{
		myarray[element - 1] = 1;
		myarray2[element - 1]++;
		myarray[element] = 0;
}


void shiftright(int myarray[], int myarray2[], int element, int size) // Shifts atom to the right
{
		myarray[element + 1] = 1;
		myarray2[element + 1]++;
		myarray[element] = 0;
}

void PressKeyToContinue()
{
	int c;
	printf("Press any key to continue... ");
	fflush(stdout);
	do c = getchar(); while ((c != '\n') && (c != EOF));
}