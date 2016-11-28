#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <random>

using namespace std;

void shiftright(int myarray[], int myarray2[], int element, int size, int boundary);
void shiftleft(int myarray[], int myarray2[], int element, int size, int boundary);

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
		int a, b, c, d, e, n, Rest;
		double plain, x, xsum, xsqsum, r, max, p;
		//// Specifying the parameters of the simulation ////
		cout << "First, enter parameters for your 1d lattice:" << endl;
		cout << "Size of the lattice: "; // The size of one dimensional array
		cin >> L;
		cout << "\nNumber of gas particles: ";
		cin >> N;
		cout << "\nSet the probability for atom to go left: ";
		cin >> p; // Sets the probability to go to the left (and to the rigth: 1 - p )
		cout << "\nChoose boundary conditions type (type '0' for periodic, '1' for reflecting b.c.): ";
		cin >> bc; // I will fix this and make user be able to type 'periodic' or 'reflecting'
		cout << "\nHow many steps for each walk: ";
		cin >> dt;
		cout << "\nHow many times you want to simulate the walk: ";
		cin >> nwalks;
		if (N >= L) {
			cout << "Invalid input! In lattice gas model number of particles cannot be greated then number of lattice cells!" << endl;
			system("Pause");
			return 1;
		}

		//// Random walk algorithm ////
		// Mersenne Twister RNG //
		random_device rd{};
		mt19937 RNG{ rd() };
		uniform_int_distribution<int> dist(0, L - 1);
		uniform_real_distribution<double> Rand{ 0.0, 1.0 };
		// rand() RNG //  -- Bad one
		// srand((unsigned)time(NULL));
		int* Lattice = new int[L]; // Our one dimensional lattice
		int* Counter = new int[L]; // Here we will store info how many times each cell were occupied during the single run
		double* dblCounter = new double[L]; // Will store here info about all runs
		double* Prob = new double[L]; // Array with probabilities
		xsum = 0; xsqsum = 0;
		for (n = 0; n < L; n++) // Setting counter to zero
		{
			dblCounter[n] = 0; // 
		}

		//// Beginning of the simulation ////
		for (iwalk = 0; iwalk < nwalks; iwalk++) { // Number of runs
												   // Setting all lattice cells to be empty first
			for (n = 0; n < L; n++)
			{
				Lattice[n] = 0;
			}
			for (n = 0; n < L; n++) // Setting counter to zero
			{
				Counter[n] = 0; // 
			}
			// Filling the lattice with particles in random fashion
			// First particle will be placed at the center of the array to compute the probability distribution vs displacement
			Lattice[L / 2] = 1;
			// Counter[L / 2] = 1;
			for (n = 1; n < N; n++) {
				//c = rand() % L; // Chooses a random cell in our 1d array.
				//c = (int)((L)* (rand() / (RAND_MAX + 1.0)));
				a = dist(RNG);
				// Check if random number generator working properly
				/*for (n = 0; n < 10000; n++) {
				a = dist(RNG);
				Counter[a]++;
				}
				cout << "\n\nSee if randomizator works properly:" << endl;
				for (n = 0; n < L; n++) {
				cout << Counter[n] << "  " << endl;
				}*/
				//Result for rand() -- Don't like the way rand() distribute particles. Will try to find a substitution for it soon.

				if (Lattice[a] != 1) {  // Check if cell is already occupied
					Lattice[a] = 1;
					Counter[a] ++;
				}
				else { // Go back if it is occupied
					n--;
				}
			}

			// Displaying the initial positions
			cout << "\nInitial positions of the attoms in the lattice:\n";
			for (n = 0; n < L; n++) {
				cout << Lattice[n];
			}cout << endl;

			// Beginning of the single run
			for (istep = 0; istep < dt; istep++) { // Number of steps for each run
				for (iter = 0; iter < N; iter++) { // Number of iterations for each time step
												   // Picking the random particle in the array
												   // c = (int)((L - 1 - 0 + 1)* (rand() / (RAND_MAX + 1.0)));  -- Bad way, but working
					c = dist(RNG);
					if (Lattice[c] != 0) { // Working only with non-empty cells
										   // Throwing the due. Generates a random number between 0 and 1. If number is less than probability p, go to the left, otherwise go to the right
										   //r = ((double)rand() / (RAND_MAX));
						r = Rand(RNG);
						if (r < p) {
							if (Lattice[c - 1] != 1 && c != 0) {
								shiftleft(Lattice, Counter, c, L, bc);
							}
							else { // If the left neighbor cell is occupied, go to the right
								if (Lattice[c - 1] != 0 && c != 0) {
									shiftright(Lattice, Counter, c, L, bc);
								}
								else { // If we have the first cell, treat separately
									if (Lattice[L - 1] != 1) {
										shiftleft(Lattice, Counter, 0, L, bc);
									}
								}
							}
						}
						else { // For r > p, we shift our particle to the right neighbor cell
							if (Lattice[c + 1] != 1 && n != L - 1) {
								shiftright(Lattice, Counter, c, L, bc);
							}
							else { // If the right neighbor cell is occupied, go to the left
								if (Lattice[c + 1] != 0 && c != L - 1) {
									shiftleft(Lattice, Counter, c, L, bc);
								}
								else { // If we have the last cell, treat separately
									shiftright(Lattice, Counter, L - 1, L, bc);
								}
							}
						}
					}
					else {
						iter--;
					}
				}
				// End of the step cycle

				// Displaying change of the positions

				for (n = 0; n < L; n++) {
					cout << Lattice[n];
				}
				cout << endl;
			}

			/*else { // For systems without exclusion
			for (istep = 0; istep < dt; istep++) { // Number of steps for each run
			// Filling the empty lattice sites with particles in random fashion
			for (n = 0; n < N; n++) {
			//c = rand() % L; // Chooses a random cell in our 1d array.
			// j = min + (int) ((max-min+1) * (rand() / (RAND_MAX + 1.0))); // Just for reference
			c = (int)((L)* (rand() / (RAND_MAX + 1.0)));
			// Check if random number generator working properly
			if (Lattice[c] == 0) {  // Check if cell is already occupied
			Lattice[c] = 1;
			Counter[c] ++;
			}
			else { // Go back if it is occupied
			n--;
			}
			}

			for (iter = 0; iter < L; iter++) { // Number of iterations for each time step
			// Picking the random particle in the array
			c = (int)((L)* (rand() / (RAND_MAX + 1.0)));
			if (Lattice[c] != 0) { // Working only with non-empty cells
			// Throwing the due. Generates a random number between 0 and 1. If number is less than probability p, go to the left, otherwise go to the right
			r = ((double)rand() / (RAND_MAX));
			if (r < p) {
			shiftleft(Lattice, Counter, c, L, bc);
			}
			}
			else { // For r > p, we shift our particle to the right neighbor cell
			shiftright(Lattice, Counter, c, L, bc);
			}
			}
			}
			// End of the step cycle

			// Displaying change of the positions
			for (n = 0; n < L; n++) {
			cout << Lattice[n];
			}
			cout << endl;
			}*/

			// End of the run

			//Computing output parameters
			// Probability distribution
			for (n = 0; n < L; n++) {
				dblCounter[n] += Counter[n];
			}

			// Displacement -- !!! Works only for a single particle in a whole lattice!!! N must equal to 1
			for (n = 0; n < L; n++) {
				if (Lattice[n] != 0) {
					x = n - L / 2;
					xsum += x;
					xsqsum += x*x; // Accumulates squared displacement
				}
			}

			Rest = 0;
			for (n = 0; n < L; n++) {
				Rest += Lattice[n];
			}
			cout << "\nTotal particles number: before " << N << ",  after  " << Rest << endl;
			if (N != Rest) {
				cout << "Partice lost during the run!!!" << endl;
			}
		}


		// End of the simulation

		//// Computing output ////
		// Probability distribution

		for (n = 0; n < L; n++) {
			Prob[n] = dblCounter[n] / (nwalks*dt*N);
		}
		plain = 1.0 / (L); // Plain probability distribution

						   //// Displaying the output ////
						   /*cout << "\nNumber of times each cell was visited during the simulation: \n";
						   for (n = 0; n < L; n++) {
						   cout << Counter[n] << "  ";
						   }*/
		if (L < 51) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
			}
		}
		if (51 < L&&L < 501) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
				n += 10;
			}
		}
		if (501 < L&&L < 5001) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
				n += 100;
			}
		}
		if (5001 < L&&L < 50000) {
			cout << "\n\nProbability distribution to find atom in each cell: \n";
			for (n = 0; n < L; n++) {
				cout << Prob[n] << scientific << endl;
				n += 1000;
			}
		}
		cout << "\nCompare with the flat distribution 1/L: " << plain << endl << endl;
		cout << "\nAverage displacement: " << xsum / nwalks << endl;
		cout << endl << "\nMean squared displacement: " << xsqsum / nwalks << endl;
		cout << "\nCompare with the number of timesteps: " << dt << endl;

		//// Exporting the output to the text file ////
		// Header
		out_stream << "This is a document that stores the output from Random Walk program. Here in particular we work with 1d lattice gas.\n\nThe simulation input parameters:\n";
		// Simulation parameters
		out_stream << "Lattice size: " << L << " cells" << "\nNumber of atoms: " << N << "\nProbability to go left: " << p << "\nBoundary conditions (0 stands for periodic, and 1 for reflecting b.c.): " << bc << "\nNumber of walks: " << nwalks << "\nNumber of steps in each walk: " << dt << "\nMean displacement: " << xsum / nwalks << "\nMean square displacement" << xsqsum / nwalks << ".\n\n Cell number | Probability to find an atom\n";
		// Output
		for (n = 0; n < L; n++) {
			out_stream << Prob[n] << "		" << n - L / 2 << endl;
		}
		//cout << "Doing random walk to " << N << " steps" << " and " << nwalks << " times, we have" << endl;
		//cout << "The average displacement " << xsum / nwalks << "\nThe mean squared displacement " << xsqsum / nwalks << endl;
		out_stream.close();

		system("Pause");
		return 0;
	}
}


void shiftleft(int myarray[], int myarray2[], int element, int size, int boundary) // Shifts atom to the left
{
	// For periodic boundary conditions
	if (boundary == 0) {
		if (element != 0) { // If atom is NOT in the first cell
			if (myarray[element - 1] != 1) { // If left neighbor cell is not occupied
				myarray[element - 1] = 1;
				myarray2[element - 1]++;
				myarray[element] = 0;
			}
		}
		else { // If atom is in the first cell
			if (myarray[size - 1] != 1) { // If last cell is not occupied
				myarray[size - 1] = 1;
				myarray2[size - 1]++;
				myarray[element] = 0;
			}
			else { // If last cell is occupied
				if (myarray[element + 1] != 1) {
					myarray[element + 1] = 1;
					myarray2[element + 1]++;
					myarray[element] = 0;
				}
			}
		}
	}

	// For reflecting boundary conditions
	if (boundary == 1) {
		if (element != 0) { // If atom is NOT in the first cell
			if (myarray[element - 1] != 1) { // If left neighbor cell is empty
				myarray[element - 1] = 1;
				myarray2[element - 1]++;
				myarray[element] = 0;
			}
		}
		else { // Reflecting back to the right
			if (myarray[element + 1] != 1) {
				myarray[element + 1] = 1;
				myarray2[element + 1]++;
				myarray[element] = 0;
			}
		}
	}

	// For non exclusive system with reflecting bc
	if (boundary == 2) {
		if (element != 0) { // If atom is NOT in the first cell
			myarray[element - 1] = myarray[element - 1] + myarray[element];
			myarray2[element - 1] = myarray2[element - 1] + myarray[element];
			myarray[element] = 0;

		}
		else { // Reflecting back to the right
			myarray[element + 1] = myarray[element + 1] + myarray[element];
			myarray2[element + 1] = myarray2[element + 1] + myarray[element];
			myarray[element] = 0;
		}
	}

}


void shiftright(int myarray[], int myarray2[], int element, int size, int boundary) // Shifts atom to the right
{
	// For periodic boundary conditions
	if (boundary == 0) {
		if (element != size - 1) { // If atom is NOT in the last cell
			if (myarray[element + 1] != 1) { // If right neighbor cell is empty
				myarray[element + 1] = 1;
				myarray2[element + 1]++;
				myarray[element] = 0;
			}
		}
		else { // If atom is in the last cell
			if (myarray[0] != 1) { // If first cell is NOT occupied - go right
				myarray[0] = 1;
				myarray2[0]++;
				myarray[element] = 0;
			}
			else { // If it is occupied - go left
				if (myarray[element - 1] != 1) {
					myarray[element - 1] = 1;
					myarray2[element - 1]++;
					myarray[element] = 0;
				}
			}
		}
	}

	// For reflecting boundary conditions
	if (boundary == 1) {
		if (element != size - 1) { // If atom is NOT in the last cell
			if (myarray[element + 1] != 1) {
				myarray[element + 1] = 1;
				myarray2[element + 1]++;
				myarray[element] = 0;
			}
		}
		else { // Reflecting back to the left
			if (myarray[element - 1] != 1) {
				myarray[element - 1] = 1;
				myarray2[element - 1]++;
				myarray[element] = 0;
			}
		}
	}

	// For non exclusive system with reflecting bc
	if (boundary == 2) {
		if (element != size - 1) { // If atom is NOT in the last cell
			myarray[element + 1] = myarray[element + 1] + myarray[element];
			myarray2[element + 1] = myarray2[element + 1] + myarray[element];
			myarray[element] = 0;

		}
		else { // Reflecting back to the left
			myarray[element - 1] = myarray[element - 1] + myarray[element];
			myarray2[element - 1] = myarray2[element - 1] + myarray[element];
			myarray[element] = 0;
		}
	}
}
