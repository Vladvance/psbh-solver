#include "genetic_algorithm.h"
#include "sbh_utils.h"

#include <cassert>
#include <iostream>
#include <chrono>
#include <numeric>
#include <random>

#include "cxxtimer.hpp"


namespace ga {

	std::random_device rd;

	genetic_algorithm::genetic_algorithm(const std::vector<oligo>& spectrum, const size_t start_oligo_idx,
		const size_t oligo_length, const size_t sequence_length) :
		spectrum_(spectrum),
		overlap_matrix_(spectrum.size(), std::vector<size_t>(spectrum.size(), 0)),
		start_oligo_idx_(start_oligo_idx),
		oligo_length_(oligo_length),
		sequence_length_(sequence_length),
		spectrum_size_(spectrum.size()),
		best_oligo_count_(sequence_length - oligo_length + 1LL),
		old_population_(population_size_, { 0, individual_t(spectrum_size_) }),
		new_population_(population_size_, { 0, individual_t(spectrum_size_) }),
		part_sum_fitness_(population_size_)
	{
		parents_count_ = static_cast<size_t>(crossover_probability_ * population_size_);
		if (parents_count_ & 1) parents_count_++; //should be even number of parents
	}

	genetic_algorithm::genetic_algorithm(const sbh_data& problem) :
		genetic_algorithm(problem.spectrum, problem.start_oligo_idx, problem.oligo_length, problem.sequence_length)
	{}

	genetic_algorithm::genetic_algorithm(const sbh_data& problem, const cxxproperties::Properties& properties) :
		genetic_algorithm(problem)
	{
		seed = properties.get<int>("seed");
		population_size_ = properties.get<size_t>("population-size");
		crossover_probability_ = properties.get<float>("crossover-probability");
		mutation_probability_ = properties.get<float>("mutation-probability");
	}

	bool genetic_algorithm::is_perfect_overlap(const uint32_t lhs, const uint32_t rhs) const
	{
		const uint32_t mask = (1 << 2 * (oligo_length_ - 1)) - 1;
		return (lhs & mask) == rhs >> 2;
	}

	void genetic_algorithm::run()
	{
		//measure time
		cxxtimer::Timer timer;
		timer.start();

		calc_overlap_matrix();

		generate_population();

		//evaluate objective function for every individual
		for (auto& ind : old_population_) {
			ind.first = objective_function(ind.second);
		}

		size_t best_value = 0;
		size_t iters_without_improvement = 0;
		size_t gen = 0;


		//main loop
		do {
			gen++;

			generation();

			//gather information about population
			struct stats stats {};
			statistics(stats);
			//report(stats, gen);
			// summary(stats, gen);
			const auto new_best_value = stats.max;

			if (new_best_value > best_value) {
				best_value = new_best_value;
				iters_without_improvement = 0;
				if (best_value == best_oligo_count_) break;
			}
			else
			{
				iters_without_improvement++;
			}

			//swap old and new populations
			old_population_.swap(new_population_);
			//} while (gen < max_generation);
		} while (gen < max_generation && iters_without_improvement < 50);

		timer.stop();

		const auto best_individual = std::max_element(old_population_.begin(), old_population_.end(), [](const auto& lhs, const auto& rhs) {return lhs.first < rhs.first; });
		// std::cout << format_individual(best_individual->second) << std::endl;
		// std::cout << "Time elapsed: " << elapsed.count() << " milliseconds" << std::endl;
		// std::cout << "Generations count: " << gen << std::endl;
		// std::cout << "Subsequent oligos count: " << best_individual->first << std::endl;
		// std::cout << "Required oligos count: " << best_oligo_count_ << std::endl;
		// std::cout << get_sequence(best_individual->second) << std::endl;
		const auto required_count = sequence_length_ - oligo_length_ + 1;
		std::printf("%s %d %lld\n", (best_individual->first == required_count) ? "Optimal" : "Feasible", best_individual->first, timer.count<std::chrono::milliseconds>());
	}


	inline void genetic_algorithm::validate(const individual_t& ind) noexcept
	{
		size_t count = 0, p = 0;
		do {
			count++;
			p = ind[p];
			assert(count <= ind.size());
		} while (p != 0);
		assert(count == ind.size());
	}

	inline std::string genetic_algorithm::get_sequence(const individual_t& ind) const
	{
		auto result(decode_n_last(spectrum_[start_oligo_idx_].seq, oligo_length_));
		auto p = start_oligo_idx_;
		auto current_sequence_length = oligo_length_;
		do {
			const auto overlap = overlap_matrix_[p][ind[p]];
			current_sequence_length += oligo_length_ - overlap;
			if (current_sequence_length > sequence_length_) break;
			//if (overlap == 0) break;
			result.append(decode_n_last(spectrum_[ind[p]].seq, oligo_length_ - overlap));
			p = ind[p];
		} while (ind[p] != start_oligo_idx_);
		return result;
	}

	inline std::string genetic_algorithm::format_individual(const individual_t& ind) const
	{
		std::string result;
		auto p = start_oligo_idx_;
		auto sequence_length = oligo_length_;
		bool is_ended = false;
		do {
			const auto overlap = overlap_matrix_[p][ind[p]];
			sequence_length += oligo_length_ - overlap;
			if (sequence_length > sequence_length_ && !is_ended) {
				is_ended = true;
				result.append("--------------------------------------------\n");
			}
			if (overlap != oligo_length_ - 1)
				result.append("Bad overlap: " + std::to_string(p) + ":" + decode_n_last(spectrum_[p].seq, oligo_length_) + "-" + std::to_string(ind[p]) + ":" + decode_n_last(spectrum_[ind[p]].seq, oligo_length_) + " = " + std::to_string(overlap) + "\n");
			//result.append(std::string(sequence_length - overlap, ' '));
			//result.append(decode_n_last(spectrum_[p].seq, oligo_length_) + "\n");
			p = ind[p];
		} while (p != start_oligo_idx_);
		return result;
	}

	//inline int genetic_algorithm::objective_function(const individual_t& ind) const
	//{
	//	std::vector <size_t> overlaps(spectrum_size_);

	//	auto p = start_oligo_idx_;
	//	size_t overlaps_idx = 0;
	//	do {
	//		//overlaps[overlaps_idx++] = calc_overlap(spectrum[p].seq, spectrum[ind[p]].seq, oligo_length);
	//		overlaps[overlaps_idx++] = overlap_matrix_[p][ind[p]];
	//		p = ind[p];
	//	} while (p != start_oligo_idx_);

	//	size_t current_end;
	//	size_t current_start = current_end = 0;
	//	size_t current_length = oligo_length_;
	//	size_t max_count = 1;

	//	for (size_t i = 0; i < overlaps.size(); ++i) {
	//		if (overlaps[i] == 0) {
	//			current_start = ++current_end;
	//			if (current_start < overlaps.size()) current_length = oligo_length_;
	//			continue;
	//		}
	//		current_length += (oligo_length_ - overlaps[i]);
	//		current_end++;
	//		while (current_length > sequence_length_)
	//			current_length -= (oligo_length_ - overlaps[current_start++]);
	//		max_count = std::max(max_count, current_end - current_start + 1);
	//	}
	//	return max_count;
	//}

	inline int genetic_algorithm::objective_function(const individual_t& ind) const noexcept
	{
		auto current_start = start_oligo_idx_;
		size_t current_end = current_start;
		size_t current_length = oligo_length_;
		size_t max_count = 1;
		size_t current_count = 1;
		bool is_end_passed = false;

		do {
			if (overlap_matrix_[current_end][ind[current_end]] == 0)
			{
				if (current_start != current_end) {
					current_length = oligo_length_;
					current_count = 1;
				}

				current_start = current_end = ind[current_end];
				continue;
			}
			current_length += oligo_length_ - overlap_matrix_[current_end][ind[current_end]];
			current_end = ind[current_end];
			current_count++;

			while (current_length > sequence_length_)
			{
				current_length -= oligo_length_ - overlap_matrix_[current_start][ind[current_start]];
				current_start = ind[current_start];
				current_count--;
			}
			max_count = std::max(max_count, current_count);
		} while (ind[current_end] != start_oligo_idx_);

		return max_count;
	}

	inline float genetic_algorithm::fitness_function(const int objective) const noexcept
	{
		return static_cast<float>(objective) / static_cast<float>(best_oligo_count_) * scaling_factor_;
	}

	inline void genetic_algorithm::calc_overlap_matrix() noexcept
	{
		for (size_t i = 0; i < spectrum_size_; ++i)
		{
			for (size_t j = 0; j < spectrum_size_; ++j)
			{
				if (i == j) continue;
				overlap_matrix_[i][j] = calc_overlap(spectrum_[i].seq, spectrum_[j].seq, oligo_length_);
			}
		}
	}

	
	// Generate temporal array of indices,
	// shuffle it and map to individual using adjacency representation
	inline void genetic_algorithm::generate_individual(individual_t& ind) noexcept
	{
		individual_t tmp(ind.size());
		std::iota(tmp.begin(), tmp.end(), 0);
		std::shuffle(tmp.begin(), tmp.end(), rd);

		auto index = ind[tmp.back()] = tmp.front();
		for (auto i = 1u; i < tmp.size(); ++i) {
			index = ind[index] = tmp[i];
		}
	}

	// Generate initial random population of individuals
	inline void genetic_algorithm::generate_population() noexcept
	{
		for (auto& ind : old_population_) {
			generate_individual(ind.second);
			ind.first = objective_function(ind.second);
		}
	}

	/**
	 * \brief Find indices of predecessors of given oligonucleotide in both parents
	 * \param parent1 first parent
	 * \param parent2 second parent
	 * \param oligo oligonucleotide which predecessors must be found
	 * \param result array of predecessors indices
	 */
	inline void genetic_algorithm::find_predecessors_in_parents(const individual_t& parent1, const individual_t& parent2,
		const uint32_t oligo, size_t result[2]) noexcept
	{
		const auto pred1 = std::find(parent1.begin(), parent1.end(), oligo);
		const auto pred2 = std::find(parent2.begin(), parent2.end(), oligo);
		result[0] = std::distance(parent1.begin(), pred1);
		result[1] = std::distance(parent2.begin(), pred2);
	}

	 // Apply mutation on given individual with probability mutation_probability_
	 // individual individual on which mutation should be applied
	inline void genetic_algorithm::mutation(individual_t& individual) const noexcept
	{
		const std::bernoulli_distribution is_mutation(mutation_probability_);

		// Return if mutation should not be applied
		if (!is_mutation(rd)) {
			return;
		}

		// pred - predecessor, succ - successor
		size_t overlap_with_pred, min_total_overlap = INT_MAX;
		auto idx_min_total_overlap = -1;
		const auto first_overlap_with_succ = overlap_with_pred = overlap_matrix_[0][individual[0]];
		auto is_overlap_with_pred_less = true;

		auto p = start_oligo_idx_;
		do {
			p = individual[p];
			const auto overlap_with_succ = calc_overlap(spectrum_[p].seq, spectrum_[individual[p]].seq, oligo_length_);
			if (overlap_with_pred + overlap_with_succ < min_total_overlap) {
				min_total_overlap = overlap_with_pred + overlap_with_succ;
				is_overlap_with_pred_less = overlap_with_pred < overlap_with_succ;
				idx_min_total_overlap = p;
			}
			overlap_with_pred = overlap_with_succ;
		} while (p != start_oligo_idx_);

		if (overlap_with_pred + first_overlap_with_succ < min_total_overlap)
			idx_min_total_overlap = 0;

		if (is_overlap_with_pred_less) {
			// Swap with predecessor
			const auto it_pred = std::find(individual.begin(), individual.end(), idx_min_total_overlap);
			const auto idx_pred = std::distance(individual.begin(), it_pred);
			const auto it_pred_pred = std::find(individual.begin(), individual.end(), idx_pred);
			const auto idx_pred_pred = std::distance(individual.begin(), it_pred_pred);
			individual[idx_pred_pred] = idx_min_total_overlap;
			individual[idx_pred] = individual[idx_min_total_overlap];
			individual[idx_min_total_overlap] = idx_pred;
		}
		else {
			// Swap with successor
			const auto it_pred = std::find(individual.begin(), individual.end(), idx_min_total_overlap);
			const auto idx_pred = std::distance(individual.begin(), it_pred);
			individual[idx_pred] = individual[idx_min_total_overlap];
			const auto tmp = individual[idx_min_total_overlap];
			individual[idx_min_total_overlap] = individual[tmp];
			individual[tmp] = idx_min_total_overlap;
		}
	}

	/**
	 * \brief Apply crossover on two selected individuals with probability crossover_probability_.
	 * Produces one offspring and applies mutation on it.
	 * \param parent1 first individual participating in crossover
	 * \param parent2 second individual participating in crossover
	 * \param offspring offspring produced after crossover. If crossover is not applied, random parent becomes offspring
	 */
	inline void genetic_algorithm::crossover(const individual_t& parent1, const individual_t& parent2, individual_t& offspring) const
	{
		const std::bernoulli_distribution is_crossover(crossover_probability_);

		if (!is_crossover(rd)) {

			// Randomly choose one parent to become offspring
			const std::bernoulli_distribution is_first_parent(0.5);
			if (is_first_parent(rd)) {
				offspring = parent1;
			}
			else {
				offspring = parent2;
			}
			mutation(offspring);
			return;
		}

		// Begin from random oligo
		const std::uniform_int_distribution<int> rand_idx(0, spectrum_size_ - 1);
		size_t oligo_first_idx = rand_idx(rd);
		size_t oligo_last_idx = oligo_first_idx;


		std::vector <bool> is_in_offspring(spectrum_size_, false);
		// Mark this oligo so we know it's already in offspring
		is_in_offspring[oligo_first_idx] = true;

		// Impossible index value to indicate max index hasn't been chosen yet
		const auto OVERLAP_UNDEFINED = spectrum_size_;

		for (size_t oligo_num = 0; oligo_num < spectrum_size_ - 1; ++oligo_num)
		{
			int p1_overlap = -1;
			int p2_overlap = -1;
			int max_overlap_idx = -1;

			if(!is_in_offspring[parent1[oligo_last_idx]])
				p1_overlap = overlap_matrix_[oligo_last_idx][parent1[oligo_last_idx]];
			if(!is_in_offspring[parent2[oligo_last_idx]])
				p2_overlap = overlap_matrix_[oligo_last_idx][parent2[oligo_last_idx]];
			if(p1_overlap != -1 || p2_overlap != -1) {
				if(p1_overlap > p2_overlap) {
					max_overlap_idx = parent1[oligo_last_idx];
				} else {
					max_overlap_idx = parent2[oligo_last_idx];
				}
			}

			// If both successors are in offspring, choose from spectrum
			if(max_overlap_idx == -1)
			{
				size_t max_overlap = OVERLAP_UNDEFINED;
				for(size_t i = 0; i < spectrum_size_; ++i)
				{
					if(is_in_offspring[i] || i == oligo_last_idx)
						continue;
					if(overlap_matrix_[oligo_last_idx][i] > max_overlap || max_overlap == OVERLAP_UNDEFINED)
					{
						max_overlap_idx = i;
						max_overlap = overlap_matrix_[oligo_last_idx][i];
					}
						
				}
			}
			offspring[oligo_last_idx] = max_overlap_idx;
			is_in_offspring[max_overlap_idx] = true;
			oligo_last_idx = max_overlap_idx;
		}
		offspring[oligo_last_idx] = oligo_first_idx;
		mutation(offspring);
	}

	inline size_t genetic_algorithm::partsum_select() const
	{
		const std::uniform_int_distribution<size_t> rand_sum(0, part_sum_fitness_.back());
		const size_t rand = rand_sum(rd);
		auto it = std::upper_bound(part_sum_fitness_.begin(), part_sum_fitness_.end(), rand);
		if(it == part_sum_fitness_.end())
			--it;
		return std::distance(part_sum_fitness_.begin(), it);
	}

	inline void genetic_algorithm::calc_part_sum_fitness()
	{
		part_sum_fitness_[0] = old_population_[0].first;
		for (size_t i = 1; i < population_size_; ++i) {
			part_sum_fitness_[i] = part_sum_fitness_[i - 1] + old_population_[i].first;
		}
	}

	inline void genetic_algorithm::generation() {

		calc_part_sum_fitness();

		// Keep best individual
		const auto max = max_element(old_population_.begin(), old_population_.end(), [](const auto& lhs, const auto& rhs)
			{
				return lhs.first < rhs.first;
			});
		new_population_[0] = *max;

		for (auto it = std::next(new_population_.begin()); it != new_population_.end(); ++it){
			const auto mate1_idx = partsum_select();
			const auto mate2_idx = partsum_select();

			crossover(old_population_[mate1_idx].second, old_population_[mate2_idx].second, it->second);
			it->first = objective_function(it->second);
		}
	}

	inline void genetic_algorithm::statistics(stats& stats) const
	{
		const auto minmax_pair = minmax_element(old_population_.begin(), old_population_.end(), [](const std::pair<int, std::vector<size_t>>& lhs, const std::pair<int, individual_t>& rhs) {return lhs.first < rhs.first; });
		const long sum = accumulate(old_population_.begin(), old_population_.end(), 0.0f, [](const long curr_sum, const std::pair<int, individual_t>& rhs) {return curr_sum + rhs.first; });
		stats.min = minmax_pair.first->first;
		stats.max = minmax_pair.second->first;
		stats.max_idx = std::distance(old_population_.begin(), minmax_pair.second);
		stats.avg = static_cast<float>(sum) / static_cast<float>(population_size_);
	}

	inline void genetic_algorithm::report(const stats& stats, const int gen) const
	{
		const auto ind_string_length = format_individual(old_population_[0].second).size();
		std::printf("------------------------------------------------------------------------\n");
		std::printf("	                    Population report                                   \n");
		std::printf("                              Generation %3d                            \n", gen);
		std::printf("------------------------------------------------------------------------\n");
		std::printf("   #               individual                        objective  fitness \n");
		std::printf("------------------------------------------------------------------------\n");
		for (size_t i = 0; i < population_size_; ++i) {
			const int objective = old_population_[i].first;
			// may add individual display
			std::printf(" %4llu %45s     %3d       %3.2f\n", i, "Individual", objective, fitness_function(objective));
		}
	}

	inline void genetic_algorithm::summary(const stats& stats, const int gen)
	{
		std::printf("------------------------------------------------------------------------\n");
		std::printf("Generation no: %2d,  max = %3llu    min = %3llu  avg = %3.3f               \n", gen, stats.max, stats.min, stats.avg);
		std::printf("------------------------------------------------------------------------\n\n");
	}
}

