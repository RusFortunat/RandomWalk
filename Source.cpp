#include <iostream>
#include <iomanip>
#include <random>

int main()
{
    std::cout << "Welcome to random walk simulation program" << std::endl << std::endl;

    // ==================================================
    // Simulation parameters input
    // ==================================================

    int latticeSize, stepsQuantity, repeatsQuantity;
    double rightProbability;

    std::cout << "Please, enter simulation parameters" << std::endl;

    std::cout << "Size of lattice: "; std::cin >> latticeSize;
    std::cout << "Probability for particle to go right: "; std::cin >> rightProbability;
    std::cout << "Simulation steps quantity: "; std::cin >> stepsQuantity;
    std::cout << "Simulation repeats quantity: "; std::cin >> repeatsQuantity;

    std::cout << std::endl;

    // ==================================================
    // Simulation parameters check
    // ==================================================

    auto errorInParameters = false;

    if (latticeSize < 1)
    {
        std::cout << "Error in simulation parameters: size of lattice must be greather or equal than 1" << std::endl;
        errorInParameters = true;
    }
    if (rightProbability <= 0 || rightProbability >= 1.0)
    {
        std::cout << "Error in simulation parameters: probability for particle to go right must be between 0 (excluding) and 1 (excluding)" << std::endl;
        errorInParameters = true;
    }
    if (stepsQuantity < 1 || stepsQuantity >= latticeSize / 2)
    {
        std::cout << "Error in simulation parameters: simulation steps quantity must be between 1 (including) and half of lattice size (excluding)" << std::endl;
        errorInParameters = true;
    }
    if (repeatsQuantity < 1)
    {
        std::cout << "Error in simulation parameters: simulation repeats quantity must be greather or equal than 1" << std::endl;
        errorInParameters = true;
    }

    if (errorInParameters)
    {
        system("PAUSE");
        return 0;
    }

    // ==================================================
    // Random numbers generator creating
    // ==================================================

    std::random_device randomDevice {};
    std::mt19937 randomGenerator { randomDevice() };

    std::uniform_real_distribution<double> numbersDistribution { 0.0, 1.0 };

    // ==================================================
    // Memory allocation
    // ==================================================

    auto probability = new double[latticeSize];
    for (auto i = 0; i < latticeSize; i++)
    {
        probability[i] = 0.0;
    }

    // ==================================================
    // Simulation
    // ==================================================

    std::cout << "Simulation started..." << std::endl;

    for (auto repeatNumber = 0; repeatNumber < repeatsQuantity; repeatNumber++)
    {
        auto particlePosition = latticeSize / 2;
        probability[latticeSize / 2]++;

        for (auto stepNumber = 0; stepNumber < stepsQuantity; stepNumber++)
        {
            auto randomNumber = numbersDistribution(randomGenerator);
            if (randomNumber <= rightProbability)
            {
                particlePosition++;
            }
            else
            {
                particlePosition--;
            }

            probability[particlePosition]++;
        }

        std::cout << std::fixed << std::setprecision(2) << std::setw(7) << static_cast<double>(repeatNumber + 1) / repeatsQuantity * 100 << "% competed" << std::endl;
    }

    std::cout << "Simulation finished" << std::endl << std::endl;

    // ==================================================
    // Simulation results calculating
    // ==================================================

    for (auto i = 0; i < latticeSize; i++)
    {
        probability[i] = probability[i] / stepsQuantity / repeatsQuantity;
    }

    // ==================================================
    // Results parameters input
    // ==================================================

    int resultsStep;

    std::cout << "Please, enter results parameters" << std::endl;
    std::cout << "Step in lattice to output results: "; std::cin >> resultsStep;
    std::cout << std::endl;

    // ==================================================
    // Results parameters check
    // ==================================================

    if (resultsStep < 1 || resultsStep >= latticeSize / 2)
    {
        std::cout << "Error in results parameters: step in lattice to output results must be between 1 (including) and half of lattice size (excluding)" << std::endl;
        system("PAUSE");
        return 0;
    }

    // ==================================================
    // Simulation results output
    // ==================================================

    std::cout << "Propability of finding a particle in a cells of lattice" << std::endl;
    for (auto i = 0; i < latticeSize; i += resultsStep)
    {
        std::cout << std::fixed << std::setprecision(2) << std::setw(7) << probability[i] * 100 << "%" << std::endl;
    }
    std::cout << std::endl;

    // ==================================================
    // Memory deallocation
    // ==================================================

    delete[] probability;
    
    system("PAUSE");
    return 0;
}