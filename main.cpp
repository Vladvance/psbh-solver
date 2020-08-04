#include <ctime>
#include <iostream>
#include <chrono>

#include "sbh_utils.h"
#include "genetic_algorithm.h"
#include "exact_algorithm.h"


int main()
{
	srand(time(nullptr));

	//const struct sbh_data problem = load_problem_hackerrank(std::cin);
	//std::ifstream is("./in.txt");
	//if(is.fail())
	//{
	//	std::cout << "No input file provided." << std::endl;
	//	exit(EXIT_FAILURE);
	//}
	const auto problem = load_problem_xml(std::cin);

	ga::genetic_algorithm ga(problem);
	//ea::exact_algorithm ea(problem);

	ga.run();

	return 0;
}

