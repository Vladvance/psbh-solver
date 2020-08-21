#pragma once

#include <string>
#include <vector>


#include "cxxproperties.hpp"
#include "sbh_utils.h"

namespace ga {

	constexpr int max_generation = 5000;

	struct stats {
		size_t min, max;
		float avg;
		size_t max_idx;
	};

	typedef std::vector<size_t> individual_t;
	typedef std::vector<std::pair<int, individual_t>> population_t;

	class genetic_algorithm {

	public:
		genetic_algorithm(const std::vector<oligo>& spectrum, size_t start_oligo_idx, size_t oligo_length,
		                  size_t sequence_length);
		genetic_algorithm(const sbh_data& problem);
		genetic_algorithm(const sbh_data& problem, const cxxproperties::Properties& properties);
		bool is_perfect_overlap(uint32_t lhs, uint32_t rhs) const;
		void preprocess();
		void run();

	private:
		const std::vector<oligo> spectrum_;
		std::vector< std::vector <size_t> > overlap_matrix_;
		const size_t start_oligo_idx_;
		const size_t oligo_length_;
		const size_t sequence_length_;

		const size_t spectrum_size_;
		const size_t best_oligo_count_;
		size_t parents_count_;

		int seed = -1;
		size_t population_size_ = 50;
		float mutation_probability_ = 0.3f;
		float crossover_probability_ = .9f;
		const float scaling_factor_ = 1.5f;

		population_t old_population_;
		population_t new_population_;
		std::vector<int> part_sum_fitness_;

		void calc_overlap_matrix() noexcept;
		static void generate_individual(individual_t& ind) noexcept;

		void generate_population() noexcept;

		static void find_predecessors_in_parents(const individual_t& parent1, const individual_t& parent2, uint32_t oligo, size_t result[2]) noexcept;

		static void validate(const individual_t& ind) noexcept;
		std::string get_sequence(const individual_t& ind) const ;
		std::string format_individual(const individual_t& ind) const ;
		int objective_function(const individual_t& ind) const noexcept;
		float fitness_function(int objective) const noexcept;
		void mutation(individual_t& individual) const noexcept;
		void crossover(const individual_t& parent1, const individual_t& parent2, individual_t& offspring) const ;

		//stochastic remainder sampling with replacement
		//std::vector<int> preselect(float avg_fitness) const;
		//static int select(std::vector<int>& choices, size_t& nremain);
		size_t partsum_select() const;
		void calc_part_sum_fitness();
		void generation();
		void statistics(stats& stats) const;
		void report(const stats& stats, int gen) const;
		static void summary(const stats& stats, int gen);
	};
}

