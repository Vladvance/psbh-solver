#include <list>

#include "exact_algorithm.h"


#include <iostream>
#include <numeric>

#include "cxxtimer.hpp"

namespace ea
{
	struct comp_seq
	{
		bool operator() (const oligo& lhs, const uint32_t seq) const { return (lhs.seq >> 2) < seq; }
		bool operator() (const uint32_t seq, const oligo& rhs) const { return seq < (rhs.seq >> 2); }
	};

	void exact_algorithm::generate_graph()
	{
		for (int i = 0; i < spectrum_size_; ++i)
		{
			for (int j = 0; j < spectrum_size_; ++j)
			{
				if(i == j) continue;
				if(calc_overlap(spectrum_[i].seq, spectrum_[j].seq, oligo_length_) == oligo_length_ - 1)
					graph_[i].push_back(j);
			}
		}
	}

	// void exact_algorithm::generate_graph() {
	// 	const uint32_t mask = (1 << 2 * (oligo_length_ - 1)) - 1;
	// 	// For each oligo in spectrum find all its successors using binary search
	// 	for (size_t i = 0; i < spectrum_size_; ++i)
	// 	{
	// 		const auto seq = spectrum_[i].seq & mask;
	// 		auto it = std::lower_bound(spectrum_.begin(), spectrum_.end(), spectrum_[i].seq & mask, comp_seq{});
	//
	// 		// Skip if oligo doesn't have successors 
	// 		if (it == spectrum_.end()) continue;
	//
	// 		// After first successor found, add other successors if present
	// 		while (it != spectrum_.end() && (spectrum_[i].seq & mask) == (it->seq >> 2))
	// 		{
	// 			const int j = std::distance(spectrum_.begin(), it);
	// 				graph_[i].push_back(j);
	// 			// if (spectrum_[i].posL < spectrum_[j].posH && spectrum_[j].posL < spectrum_[i].posH)
	// 			++it;
	// 		}
	// 	}
	// }

	// Struct holds adjacent representation of solution
	typedef struct solution
	{
		solution(const size_t last_idx, const size_t oligo_count, const size_t spectrum_size, const size_t default_value)
			: last_idx(last_idx),
			  oligo_count(oligo_count),
			  data(spectrum_size, default_value)
		{ }

		size_t last_idx;
		size_t oligo_count;
		std::vector<size_t> data;
	} solution_t;

	void exact_algorithm::run() {
		cxxtimer::Timer timer;
		timer.start();

		generate_graph();

		const size_t INDEX_UNDEFINED = spectrum_size_;
		const size_t required_count = sequence_length_ - oligo_length_ + 1;

		std::vector <solution_t> solutions(1, solution_t(start_oligo_idx_, 1, spectrum_size_, INDEX_UNDEFINED));
		size_t sol_idx = 0;
		while(sol_idx < solutions.size()) {
			// Begin construction from last added oligo
			size_t idx = solutions[sol_idx].last_idx;
			
			// If last oligo in solution has child in graph
			// AND this child doesn't introduce subcycle
			// AND required oligo count haven't been reached yet
			// - append next oligo to solution
			while (!graph_[idx].empty() &&
					solutions[sol_idx].data[idx] == INDEX_UNDEFINED &&
					solutions[sol_idx].oligo_count <= required_count) 
			{
				if (graph_[idx].size() > 1)
				{
					// If there is more then one child for last oligo, for each subsequent child
					// create copy of current solution and append that child to copy
					for (auto it = std::next(graph_[idx].begin()); it != graph_[idx].end(); ++it)
					{
						if(is_position_used && (solutions[sol_idx].oligo_count + 1 < spectrum_[*it].posL ||
						   spectrum_[*it].posH < solutions[sol_idx].oligo_count + 1)) 
							continue;
						solutions.emplace_back(solutions[sol_idx]);
						solutions.back().data[idx] = *it;
						solutions.back().last_idx = *it;
						solutions.back().oligo_count++;
					}
				}
				const size_t next_oligo_idx = graph_[idx].front();
				// If next_oligo doesn't fit, stop construction and move to next solution
				if (is_position_used && (solutions[sol_idx].oligo_count + 1 < spectrum_[next_oligo_idx].posL ||
					spectrum_[next_oligo_idx].posH < solutions[sol_idx].oligo_count + 1))
					break;;
				solutions[sol_idx].data[idx] = next_oligo_idx;
				solutions[sol_idx].last_idx = next_oligo_idx;
				solutions[sol_idx].oligo_count++;
				idx = next_oligo_idx;
			}
			if(solutions[sol_idx].oligo_count == required_count) break;
			sol_idx++;
		}
		if(sol_idx == solutions.size()) sol_idx--;

		// Reconstruct result DNA string from solution
		std::string result(decode_n_last(spectrum_[start_oligo_idx_].seq, oligo_length_));
		result.reserve(required_count);
		size_t p = start_oligo_idx_;
		size_t counter = 1;
		// std::printf("Pos = %4llu; Idx = %4llu posL = %4llu; posH = %4llu; %s; nbranch = %llu%c\n", counter, p, spectrum_[p].posL, spectrum_[p].posH, decode_n_last(spectrum_[p].seq, oligo_length_).c_str(), graph_[p].size(), graph_[p].size()>1 ? '+' : '-');
		do {
			counter++;
			p = solutions[sol_idx].data[p];
			// std::printf("Pos = %4llu; Idx = %4llu posL = %4llu; posH = %4llu; %s; nbranch = %llu%c\n", counter, p, spectrum_[p].posL, spectrum_[p].posH, decode_n_last(spectrum_[p].seq, oligo_length_).c_str(), graph_[p].size(), graph_[p].size()>1 ? '+' : '-');
			result.append(decode_n_last(spectrum_[p].seq, 1));
		} while (p != solutions[sol_idx].last_idx);

		timer.stop();
		std::printf("%s %llu %lld\n", (solutions[sol_idx].oligo_count == required_count) ? "Optimal" : "Feasible", solutions[sol_idx].oligo_count, timer.count<std::chrono::milliseconds>());
		// std::printf("Time: %lld milliseconds.\n", timer.count<std::chrono::milliseconds>());
		// std::printf("Final sequence:\n");
		// std::cout << result;
		return;
	}
}
