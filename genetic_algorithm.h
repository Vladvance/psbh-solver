#pragma once

#include <string>
#include <vector>

#include "sbh_utils.h"

namespace ga {

	constexpr int max_generation = 5000;

	struct stats {
		size_t min, max;
		float avg;
		size_t max_idx;
	};


	typedef std::vector<unsigned> individual_t;
	typedef std::vector<std::pair<int, individual_t>> population_t;


	class genetic_algorithm {

	public:
		genetic_algorithm(const std::vector<oligo>& spectrum, unsigned start_oligo_idx, int oligo_length,
		                  int sequence_length);
		genetic_algorithm(const sbh_data& problem);
		void run();

	private:
		const std::vector<oligo> spectrum_;
		std::vector< std::vector <size_t> > overlap_matrix_;
		const size_t start_oligo_idx_;
		const size_t oligo_length_;
		const size_t sequence_length_;

		const int spectrum_size_;
		const int best_oligo_count_;
		const size_t population_size_;
		size_t parents_count_;

		const float mutation_probability_ = 0.001f;
		const float crossover_probability_ = 0.9f;
		const float scaling_factor_ = 1.5f;


		population_t old_population_;
		population_t new_population_;
		std::vector<float> part_sum_fitness_;

		void calc_overlap_matrix();

		static void generate_individual(individual_t& ind);
		void generate_population();

		static void find_predecessors_in_parents(const individual_t& parent1, const individual_t& parent2, uint32_t oligo, unsigned result[2]);

		static void validate(const individual_t& ind);
		std::string get_sequence(const individual_t& ind) const;
		std::string format_individual(const individual_t& ind) const;
		int objective_function(const individual_t& ind) const;
		float fitness_function(int objective) const;
		void mutation(individual_t& individual) const;
		void crossover(const individual_t& parent1, const individual_t& parent2, individual_t& offspring) const;

		//stochastic remainder sampling with replacement
		//std::vector<int> preselect(float avg_fitness) const;
		//static int select(std::vector<int>& choices, size_t& nremain);
		int partsum_select() const;

		void calc_part_sum_fitness(float avg_fitness);
		void generation();
		void statistics(stats& stats) const;
		void report(stats& stats, int gen) const;

	};
}

