#include <list>

#include "exact_algorithm.h"
#include "cxxtimer.hpp"

namespace ea
{
	struct comp_seq
	{
		bool operator() (const oligo& lhs, const uint32_t seq) const { return (lhs.seq >> 2) < seq; }
		bool operator() (const uint32_t seq, const oligo& rhs) const { return seq < (rhs.seq >> 2); }
	};

	void exact_algorithm::generate_complete_graph() {
		const uint32_t mask = (1 << 2 * (oligo_length_ - 1)) - 1;
		//std::vector <std::vector<size_t> > graph(spectrum_size_);
		for (size_t i = 0; i < spectrum_size_; ++i)
		{
			const auto seq = spectrum_[i].seq & mask;
			auto it = std::lower_bound(spectrum_.begin(), spectrum_.end(), seq, comp_seq{});
			if (it == spectrum_.end()) continue;
			while (it != spectrum_.end() && seq == (it->seq >> 2))
			{
				const int j = std::distance(spectrum_.begin(), it);
				if (spectrum_[i].posL < spectrum_[j].posH && spectrum_[i].posH > spectrum_[j].posL)
					graph_[i].push_back(j);
				++it;
			}
		}
	}

	void exact_algorithm::run() {
		cxxtimer::Timer timer;

		timer.start();
		generate_complete_graph();
		std::vector < std::vector <size_t> > solutions;
		solutions.push_back({ start_oligo_idx_ });
		for (size_t i = 1; i < iterations_; ++i) {
			std::vector < std::vector <size_t> > new_solutions;
			for (const auto& solution : solutions) {
				const auto last_oligo_idx = solution.back();
				for (const auto successor_idx : graph_[last_oligo_idx])
				{
					if (spectrum_[successor_idx].posL <= i + 1 && i <= spectrum_[successor_idx].posH)
					{
						new_solutions.emplace_back(solution);
						new_solutions.back().push_back(successor_idx);
					}
				}
			}
			solutions = std::move(new_solutions);
		}
		timer.start();
		std::printf("Time: %lld milliseconds.\n", timer.count<std::chrono::milliseconds>());
		std::printf("Final sequence/s: \n");

		for (const auto& solution : solutions)
		{
			if (solution.size() == iterations_)
			{
				std::string result(decode_n_last(spectrum_[solution[0]].seq, oligo_length_));

				for (size_t j = 1; j < solution.size(); ++j)
				{
					result.append(decode_n_last(spectrum_[solution[j]].seq, 1));
				}
				std::printf("%s\n", result.c_str());
			}
		}
		//	// Adjacency matrix
		//	//std::vector < std::vector <bool> > tmp_graph(spectrum_size_, std::vector<bool>(spectrum_size_, false));
		//	// Adjacency list
		//	std::vector < std::list<size_t> > tmp_graph(spectrum_size_);
		//	std::vector <size_t> added_nodes;
		//	added_nodes.reserve(spectrum_size_);
		//	added_nodes.push_back(start_oligo_idx_);
		//	//std::printf("Spectrum:\n");
		//	//for(auto& o: spectrum_)
		//	//{
		//	//	std::cout << o.posL << " " << o.posH << " " << decode_n_last(o.seq, oligo_length_) << std::endl;
		//	//}
		//	//std::printf("Pos L queue:\n");
		//	//while(!posl_queue_.empty()) {
		//	//	std::cout << posl_queue_.top().posL << " " << posl_queue_.top().posH << std::endl;
		//	//	posl_queue_.pop();
		//	//}
		//	//std::printf("Pos H queue:\n");
		//	//while(!posh_queue_.empty()) {
		//	//	std::cout << posh_queue_.top().posL << " " << posh_queue_.top().posH << std::endl;
		//	//	posh_queue_.pop();
		//	//}

		//	for(size_t i = 1; i < iterations_; ++i)
		//	{
		//		while(posl_queue_.top().posL <= i+1)
		//		{
		//			auto current_node = posl_queue_.top().index;
		//			for(auto node : added_nodes)
		//			{
		//				if(node == current_node) continue;;
		//				std::printf("Check %s : %s : %s\n", decode_n_last(spectrum_[node].seq, oligo_length_).c_str(), decode_n_last(spectrum_[current_node].seq, oligo_length_).c_str(), (graph_[node][current_node]) ? "true" : "false");
		//				if(graph_[node][current_node])
		//				{
		//					tmp_graph[node].push_back(current_node);
		//					added_nodes.push_back(current_node);
		//				}
		//				if(graph_[current_node][node])
		//				{
		//					tmp_graph[current_node].push_back(node);
		//					added_nodes.push_back(current_node);
		//				}
		//				//tmp_graph[node][current_node] = graph_[node][current_node]; 
		//			}
		//			posl_queue_.pop();
		//		}

		//		while(posh_queue_.top().posH <= i) {
		//			tmp_graph[posh_queue_.top().index].clear();
		//			for(size_t n = 0; n < spectrum_size_; ++n)
		//			{
		//				//tmp_graph[n][posh_queue_.top().index] = false;
		//				tmp_graph[n].remove(posh_queue_.top().index);

		//			}
		//			posh_queue_.pop();
		//		}
		//		std::vector< std::vector<size_t> > new_solutions;
		//		for(auto& solution : solutions)
		//		{
		//			assert(solution.empty());
		//			const auto last_oligo = solution.back();

		//			if(tmp_graph[last_oligo].empty())
		//			{
		//				new_solutions.emplace_back(solution);
		//				continue;
		//			}
		//				
		//			for(auto successor: tmp_graph[last_oligo])
		//			{
		//				new_solutions.emplace_back(solution);
		//				new_solutions.back().push_back(successor);
		//			}
		//		}
		//		solutions = std::move(new_solutions);
		//	}
		//	return;
	}
}
