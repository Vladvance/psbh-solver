#include <ctime>
#include <iostream>
#include <chrono>

#include "sbh_utils.h"
#include "genetic_algorithm.h"
#include "exact_algorithm.h"
#include "cxxopts.hpp"
#include "cxxproperties.hpp"


void parse_cli_args(int argc, char** argv, cxxproperties::Properties& properties)
{
	cxxopts::Options options("PSBH Solver", "Command line solver for DNA sequencing by hybridization with positive errors and position information available.");

	options.add_options("General")
		("f,file", "Input xml file with oligo data", cxxopts::value<std::string>()->default_value("in.txt"))
		("a,algorithm", "Set algorithm which should be used", cxxopts::value<std::string>()->default_value("exact"))
		;
	options.add_options("Exact")
		("p,use-position-info", "Set if exact algorithm should use position info")
		;
	options.add_options("Genetic")
		("z,population-size", "Set population size", cxxopts::value<size_t>()->default_value("50"))
		("c,crossover-probability", "Set crossover probability", cxxopts::value<float>()->default_value("0.9f"))
		("m,mutation-probability", "Set mutation probability", cxxopts::value<float>()->default_value("0.001f"))
		("s,seed", "Set seed to init random generator", cxxopts::value<int>())
		("h,help", "Show this help", cxxopts::value<int>())
		;

	try {
		auto parse_result = options.parse(argc, argv);
		std::string filename, algorithm;
		filename = parse_result["file"].as<std::string>();
		algorithm = parse_result["algorithm"].as<std::string>();
		properties.add("file", filename);
		properties.add("algorithm", algorithm);
		if(parse_result.count("help"))
		{
			std::cout << options.help() << std::endl;
			exit(EXIT_SUCCESS);
		}
		if(parse_result.count("file")) {
			filename = parse_result["file"].as<std::string>();
			properties.add("file", filename);
		}
		if(parse_result.count("algorithm")) {
			algorithm = parse_result["algorithm"].as<std::string>();
			properties.add("algorithm", algorithm);
		}
		if(algorithm == "exact")
		{
			properties.add("use-position-info", parse_result["use-position-info"].as<bool>());
		} else if(algorithm == "genetic")
		{
			properties.add("seed", parse_result["seed"].as<int>());
			properties.add("population-size", parse_result["population-size"].as<size_t>());
			properties.add("crossover-probability", parse_result["crossover-probability"].as<float>());
			properties.add("mutation-probability", parse_result["mutation-probability"].as<float>());
		} else
		{
			std::cerr << "Error: Incorrect algorithm provided. Use '--algorithm <algorithm>' option to specify algorithm. See --help for details." << std::endl;
			exit(EXIT_FAILURE);
		}
	} catch (cxxopts::OptionException& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char** argv)
{
	cxxproperties::Properties properties;
	parse_cli_args(argc, argv, properties);
	
	const std::string filename = properties.get("file");
	const std::string algorithm = properties.get("algorithm");
	
	std::ifstream is(filename);
	if(is.fail()) {
		std::cerr << "Error: No input file provided. Use '--file <filename>' option to specify input file. See '--help' for details." << std::endl;
		exit(EXIT_FAILURE);
	}

	const auto problem = load_problem_xml(is);

	if(algorithm == "exact") {
		ea::exact_algorithm ea(problem, properties);
		ea.run();
	} else {
		ga::genetic_algorithm ga(problem, properties);
		ga.run();
	}

	return 0;
}
