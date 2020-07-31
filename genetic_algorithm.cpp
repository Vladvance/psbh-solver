#include "genetic_algorithm.h"
#include "sbh_utils.h"

#include <cassert>
#include <iostream>
#include <chrono>
#include <numeric>
#include <random>


namespace ga {


	std::random_device rd;
	std::mt19937 gen(rd());


	genetic_algorithm::genetic_algorithm(const std::vector<oligo>& spectrum, const unsigned start_oligo_idx,
		const int oligo_length, const int sequence_length) :
		spectrum_(spectrum),
		overlap_matrix_(spectrum.size(), std::vector< size_t > (spectrum.size(),  0)),
		start_oligo_idx_(start_oligo_idx),
		oligo_length_(oligo_length),
		sequence_length_(sequence_length),
		spectrum_size_(spectrum.size()),
		best_oligo_count_(sequence_length - oligo_length + 1),
		population_size_(std::max(sequence_length/2 & ~1, 50)),
		old_population_(population_size_, std::pair<int, std::vector <unsigned> >(0, std::vector <unsigned>(spectrum_size_))),
		new_population_(population_size_, std::pair<int, std::vector <unsigned> >(0, std::vector <unsigned>(spectrum_size_))),
		part_sum_fitness_(population_size_)
	{
		parents_count_ = static_cast<unsigned>(crossover_probability_ * population_size_);
		if (parents_count_ & 1) parents_count_++; //should be even number of parents
	}

	genetic_algorithm::genetic_algorithm(const sbh_data& problem) :
		genetic_algorithm(problem.spectrum, problem.start_oligo_idx, problem.oligo_length, problem.sequence_length)
	{}


	void genetic_algorithm::run()
	{
		//measure time
		const auto start = std::chrono::steady_clock::now();

		calc_overlap_matrix();

		generate_population();

		//evaluate objective function for every individual
		for (auto& ind : old_population_) {
			ind.first = objective_function(ind.second);
		}

		auto best_value = 0;
		auto iters_without_improvement = 0;
		auto gen = 0;


		//main loop
		do {
			gen++;

			generation();

			//gather information about population
			struct stats stats{};
			statistics(stats);
			//report(stats, gen);
			const auto new_best_value = stats.max;

			if (new_best_value > best_value) {
				best_value = new_best_value;
				iters_without_improvement = 0;
				if (best_value == best_oligo_count_) break;
			} else
			{
				iters_without_improvement++;
			}

			//swap old and new populations
			old_population_.swap(new_population_);
		//} while (gen < max_generation);
		} while (gen < max_generation && iters_without_improvement < 50);

		const auto end = std::chrono::steady_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

		const auto best_individual = std::max_element(old_population_.begin(), old_population_.end(), [](const auto& lhs, const auto& rhs) {return lhs.first < rhs.first; });
		std::cout << format_individual(best_individual->second) << std::endl;
		std::cout << get_sequence(best_individual->second) << std::endl;
		std::cout << "Subsequent oligos count: " << best_individual->first << std::endl;
		std::cout << "Required oligos count: " << best_oligo_count_ << std::endl;
		std::cout << "Time elapsed: " << elapsed.count() << " milliseconds" << std::endl;
		std::cout << "Generations count: " << gen << std::endl;
	}


	inline void genetic_algorithm::validate(const individual_t& ind) {
		unsigned count = 0, p = 0;
		do {
			count++;
			p = ind[p];
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
			if(sequence_length > sequence_length_ && !is_ended) {
				is_ended = true;
				result.append("--------------------------------------------\n");
			}
			if(overlap != oligo_length_ - 1)
				result.append("Bad overlap: " + std::to_string(p) + ":" + decode_n_last(spectrum_[p].seq, oligo_length_) + "-" + std::to_string(ind[p]) + ":" + decode_n_last(spectrum_[ind[p]].seq, oligo_length_) + " = " + std::to_string(overlap) + "\n");
			//result.append(std::string(sequence_length - overlap, ' '));
			//result.append(decode_n_last(spectrum_[p].seq, oligo_length_) + "\n");
			p = ind[p];
		} while (p != start_oligo_idx_);
		return result;
	}

	inline int genetic_algorithm::objective_function(const individual_t& ind) const
	{
		std::vector <unsigned> overlaps(spectrum_size_);

		auto p = start_oligo_idx_;
		unsigned overlaps_idx = 0;
		do {
			//overlaps[overlaps_idx++] = calc_overlap(spectrum[p].seq, spectrum[ind[p]].seq, oligo_length);
			overlaps[overlaps_idx++] = overlap_matrix_[p][ind[p]];
			p = ind[p];
		} while (p != start_oligo_idx_);

		size_t current_end;
		size_t current_start = current_end = 0;
		size_t current_length = oligo_length_;
		size_t max_count = 1;

		for (size_t i = 0; i < overlaps.size(); ++i) {
			if (overlaps[i] == 0) {
				current_start = ++current_end;
				if (current_start < overlaps.size()) current_length = oligo_length_;
				continue;
			}
			current_length += (oligo_length_ - overlaps[i]);
			current_end++;
			while (current_length > sequence_length_)
				current_length -= (oligo_length_ - overlaps[current_start++]);
			max_count = std::max(max_count, current_end - current_start + 1);
		}
		return max_count;
	}

	inline float genetic_algorithm::fitness_function(const int objective) const
	{
		return static_cast<float>(objective) / static_cast<float>(best_oligo_count_) * scaling_factor_; 
	}

	inline void genetic_algorithm::calc_overlap_matrix()
	{
		for (auto i = 0; i < spectrum_size_; ++i)
		{
			for (auto j = 0; j < i; ++j)
			{
				overlap_matrix_[i][j] = calc_overlap(spectrum_[i].seq, spectrum_[j].seq, oligo_length_);
				overlap_matrix_[j][i] = calc_overlap(spectrum_[j].seq, spectrum_[i].seq, oligo_length_);
			}
			
		}
	}

	/*
	 * Generate temporal array of indices,
	 * shuffle it and map to individual using adjacency representation
	 */
	inline void genetic_algorithm::generate_individual(individual_t& ind)
	{
		individual_t tmp(ind.size());
		std::iota(tmp.begin(), tmp.end(), 0);
		std::shuffle(tmp.begin(), tmp.end(), gen);

		auto index = ind[tmp.back()] = tmp.front();
		for (auto i = 1u; i < tmp.size(); ++i) {
			index = ind[index] = tmp[i];
		}
	}

	/**
	 * \brief Generate initial random population of individuals
	 */
	inline void genetic_algorithm::generate_population()
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
		const uint32_t oligo, unsigned result[2])
	{
		const auto pred1 = std::find(parent1.begin(), parent1.end(), oligo);
		const auto pred2 = std::find(parent2.begin(), parent2.end(), oligo);
		result[0] = std::distance(parent1.begin(), pred1);
		result[1] = std::distance(parent2.begin(), pred2);
	}

	/**
	 * \brief Apply mutation on given individual with probability mutation_probability_
	 * \param individual individual on which mutation should be applied
	 */
	inline void genetic_algorithm::mutation(individual_t& individual) const
	{
		const std::bernoulli_distribution is_mutation(mutation_probability_);

		// Return if mutation should not be applied
		if (!is_mutation(gen)) {
			return;
		}

		// pred - predecessor, succ - successor
		size_t overlap_with_pred, min_total_overlap = INT_MAX;
		auto idx_min_total_overlap = -1;
		const auto first_overlap_with_succ = overlap_with_pred = overlap_matrix_[0][individual[0]];
		auto is_overlap_with_pred_less = true;

		auto p = start_oligo_idx_;
		auto size = best_oligo_count_;
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

		if (!is_crossover(gen)) {

			// Randomly choose one parent to become offspring
			const std::bernoulli_distribution is_first_parent(0.5);
			if (is_first_parent(gen)) {
				offspring = parent1;
			}
			else {
				offspring = parent2;
			}
			mutation(offspring);
			return;
		}

		// Variables to track the beginning and the ending of constructed oligo sequence
		int oligo_first_idx, oligo_last_idx;

		// Store overlaps and corresponding indices of oligos adjacent to edge oligos in both parents
		int adjacent_overlaps[4]{ 0 }; //0,1 - predecessors; 2,3 - successors
		unsigned adjacent_indices[4];

		std::vector <bool> is_in_offspring(spectrum_size_, false);

		const std::uniform_int_distribution<int> rand_idx(0, spectrum_size_ - 1);

		// Begin constructing from random oligo
		oligo_first_idx = oligo_last_idx = rand_idx(gen);

		// Mark this oligo so we know it's already in offspring
		is_in_offspring[oligo_first_idx] = true; 

		// Track how many oligos are not added yet
		int oligos_left = spectrum_size_ - 1;

		// Flag shows which side of constructed sequence was modified last
		bool is_beginning_modified = true; 

		while (oligos_left--) {
			/**
			 * If last added oligo was at the beginning, calculate overlaps between oligos
			 * predeceasing to first oligo in both parents and first oligo in constructed sequence
			 * Also calculate on first iteration
			 */
			if (is_beginning_modified) {
				find_predecessors_in_parents(parent1, parent2, oligo_first_idx, adjacent_indices);
				for (auto i : { 0, 1 }) {
					if (!is_in_offspring[adjacent_indices[i]]) { 
						//adjacent_overlaps[i] = calc_overlap(spectrum[adjacent_indices[i]].seq, spectrum[oligo_first_idx].seq, oligo_length);
						adjacent_overlaps[i] = overlap_matrix_[adjacent_indices[i]][oligo_first_idx];
					}
				}
			}

			/**
			 * If last added oligo was at the ending, calculate overlaps between last oligo
			 * in constructed sequence and oligos succeeding to last oligo in both parents
			 * Also calculate on first iteration
			 */
			if (!is_beginning_modified || oligos_left == spectrum_size_ - 2) { 
				adjacent_indices[2] = parent1[oligo_last_idx];
				adjacent_indices[3] = parent2[oligo_last_idx];
				for (auto i : { 2, 3 }) {
					if (!is_in_offspring[adjacent_indices[i]]) {
						//adjacent_overlaps[i] = calc_overlap(spectrum[oligo_last_idx].seq, spectrum[adjacent_indices[i]].seq, oligo_length);
						adjacent_overlaps[i] = overlap_matrix_[oligo_last_idx][adjacent_indices[i]];
					}
				}
			}

			auto max_overlap = -1, max_overlap_idx = -1;
			for (auto i = 0; i < 4; i++) {
				if (!is_in_offspring[adjacent_indices[i]] && adjacent_overlaps[i] > max_overlap) {
					max_overlap = adjacent_overlaps[i];
					max_overlap_idx = adjacent_indices[i];
					is_beginning_modified = i < 2;
				}
			}

			/**
			* If all adjacent oligos are already in offspring,
			* search for best fitted oligo in whole spectrum
			*/
			if (max_overlap_idx == -1) {
				
				for (auto i = 0; i < spectrum_size_; ++i) {
					if (!is_in_offspring[i]) {
							//const int tmp_overlap_pred = calc_overlap(spectrum[i].seq, spectrum[oligo_first_idx].seq, oligo_length);
							//const int tmp_overlap_succ = calc_overlap(spectrum[oligo_last_idx].seq, spectrum[i].seq, oligo_length);
							const int tmp_overlap_pred = overlap_matrix_[i][oligo_first_idx];
							const int tmp_overlap_succ = overlap_matrix_[oligo_last_idx][i];
							int tmp_overlap;

							if (tmp_overlap_pred > tmp_overlap_succ) {
								is_beginning_modified = true;
								tmp_overlap = tmp_overlap_pred;
							} else
							{
								is_beginning_modified = false;
								tmp_overlap = tmp_overlap_succ;
							}

							if (tmp_overlap > max_overlap) {
								max_overlap = tmp_overlap;
								max_overlap_idx = i;
							}
					}
				}
			}

			if (is_beginning_modified) {
				offspring[oligo_last_idx] = max_overlap_idx;
				oligo_last_idx = max_overlap_idx;
				is_in_offspring[oligo_last_idx] = true;
			}
			else {
				offspring[max_overlap_idx] = oligo_first_idx;
				oligo_first_idx = max_overlap_idx;
				is_in_offspring[oligo_first_idx] = true;
			}
		}
		offspring[oligo_last_idx] = oligo_first_idx;
		mutation(offspring);
	}


	inline int genetic_algorithm::partsum_select() const
	{
		const std::uniform_real_distribution<float> rand_sum(0, part_sum_fitness_.back());
		const float rand = rand_sum(gen);
		const auto it_start = std::upper_bound(part_sum_fitness_.begin(), part_sum_fitness_.end(), rand);
		auto idx = std::distance(part_sum_fitness_.begin(), it_start);

		if (it_start < std::prev(part_sum_fitness_.end())) {
			auto it_end = it_start;
			while (*(it_end + 1) == *(it_start)) ++it_end;
			if (it_start != it_end) {
				const int range_size = std::distance(it_start, it_end);
				const std::uniform_int_distribution<int> rand_idx(0, range_size - 1);
				idx += rand_idx(gen);
			}
		}
		return idx;
	}


	inline void genetic_algorithm::calc_part_sum_fitness(const float avg_fitness)
	{
		part_sum_fitness_[0] = fitness_function(old_population_[0].first) / avg_fitness;
		for (size_t i = 1; i < population_size_; ++i) {
			const float fitness = fitness_function(old_population_[i].first) / avg_fitness;
			part_sum_fitness_[i] = part_sum_fitness_[i - 1] + fitness;
		}
	}

	inline void genetic_algorithm::generation() {

		const float sum_fitness = std::accumulate(old_population_.begin(), old_population_.end(), .0, [this](const float& sum, const std::pair<int, individual_t>& ind) {return sum + fitness_function(ind.first); });
		const float avg_fitness = sum_fitness / population_size_;
		calc_part_sum_fitness(avg_fitness);

		// Keep best individual
		const auto max = max_element(old_population_.begin(), old_population_.end(), [](const auto& lhs, const auto& rhs)
			{
				return lhs.first < rhs.first;
			});
		new_population_[0] = *max;

		for (auto i = 1u; i < population_size_; ++i) {
			const auto mate1_idx = partsum_select();
			const auto mate2_idx = partsum_select();

			crossover(old_population_[mate1_idx].second, old_population_[mate2_idx].second, new_population_[i].second);
			new_population_[i].first = objective_function(new_population_[i].second);
		}
	}

	inline void genetic_algorithm::statistics(stats& stats) const
	{
		const auto minmax_pair = minmax_element(old_population_.begin(), old_population_.end(), [](const std::pair<int, std::vector<unsigned>>& lhs, const std::pair<int, individual_t>& rhs) {return lhs.first < rhs.first; });
		const long sum = accumulate(old_population_.begin(), old_population_.end(), 0.0f, [](const long curr_sum, const std::pair<int, individual_t>& rhs) {return curr_sum + rhs.first; });
		stats.min = minmax_pair.first->first;
		stats.max = minmax_pair.second->first;
		stats.max_idx = std::distance(old_population_.begin(), minmax_pair.second);
		stats.avg = static_cast<float>(sum) / static_cast<float>(population_size_);
	}

	inline void genetic_algorithm::report(stats& stats, const int gen) const
	{
		const auto ind_string_length = format_individual(old_population_[0].second).size();
		std::printf("------------------------------------------------------------------------\n");
		std::printf("	                    Population report                                   \n");
		std::printf("                              Generation %3d                            \n", gen);
		std::printf("------------------------------------------------------------------------\n");
		std::printf("   #               individual                        objective  fitness \n");
		std::printf("------------------------------------------------------------------------\n");
		for (auto i = 0; i < population_size_; ++i) {
			const int objective = old_population_[i].first;
			// may add individual display
			std::printf(" %4d %45s     %3d       %3.2f\n", i, "Individual", objective, fitness_function(objective));
		}
		std::printf("------------------------------------------------------------------------\n");
		std::printf("   max = %3d    min = %3d  avg = %3.3f                                  \n", stats.max, stats.min, stats.avg);
		std::printf("------------------------------------------------------------------------\n\n");
	}
}

