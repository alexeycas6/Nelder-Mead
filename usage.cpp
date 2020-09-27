#include <iostream>

#include "NelderMead.h"

int main() {

	// Rosenbrock function
	auto function = [](std::vector<double> x) {
		double sum = 0.0;
		int dimension = 3;
		for (int i = 0; i <= dimension - 2; i++) {
			double value = 100.0 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) +
				(x[i] - 1.0) * (x[i] - 1.0);
			sum += value;
		}
		return sum;
	};

	NelderMead optimization(function, { -4.0, -4.0, -4.0 }, { 4.0, 4.0, 4.0 });

	// exit by closeness simlex points
	auto res1 = optimization.getResult(false, 1.E-5);
	std::cout << "Exit by closeness simlex points" << std::endl
	<< "point: (";
	for (const auto& elem : res1.first)
		std::cout << elem << " ";
	std::cout << "), value = " << res1.second << std::endl << std::endl;

	// exit by time	
	auto res2 = optimization.getResult(true, 3.0);
	std::cout << "Exit by time" << std::endl
	<< "point: (";
	for (const auto& elem : res2.first)
		std::cout << elem << " ";
	std::cout << "), value = " << res2.second << std::endl;

	system("pause");
	return 0;
}