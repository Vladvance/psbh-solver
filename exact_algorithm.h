#pragma once
#include <list>
#include <queue>

#include "sbh_utils.h"

namespace ea {

	struct greater_then_posh {
		bool operator()(const oligo& lhs, const oligo& rhs) const { return lhs.posH > rhs.posH; }
	};

	struct greater_then_posl {
		bool operator()(const oligo& lhs, const oligo& rhs) const { return lhs.posL > rhs.posL; }
	};
	
	class exact_algorithm
	{
	private:
		std::priority_queue<oligo, std::vector<oligo>, greater_then_posh> posh_queue_;
		std::priority_queue<oligo, std::vector<oligo>, greater_then_posl> posl_queue_;
		const greater_then_posh posH_comparator_{};
		const greater_then_posl posL_comparator_{};

		const size_t sequence_length_;
		const size_t oligo_length_;
		const size_t start_oligo_idx_;
		std::vector<oligo> spectrum_;
		const size_t spectrum_size_;
		const size_t iterations_;
		std::vector < std::list<size_t> > graph_;
		void generate_complete_graph();
	public:
		void run();

		explicit exact_algorithm(const sbh_data& problem) :
			posh_queue_(posH_comparator_, problem.spectrum),
			posl_queue_(posL_comparator_, problem.spectrum),
			sequence_length_(problem.sequence_length),
			oligo_length_(problem.oligo_length),
			start_oligo_idx_(problem.start_oligo_idx),
			spectrum_(problem.spectrum),
			spectrum_size_(problem.spectrum.size()),
			iterations_(problem.sequence_length - problem.oligo_length + 1),
			graph_(problem.spectrum.size())
		{ }
	};
}
