#pragma once

#include "cxxproperties.hpp"
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
	public:
		void run();

		explicit exact_algorithm(const sbh_data& problem) :
			sequence_length_(problem.sequence_length),
			oligo_length_(problem.oligo_length),
			start_oligo_idx_(problem.start_oligo_idx),
			spectrum_(problem.spectrum),
			spectrum_size_(problem.spectrum.size()),
			graph_(problem.spectrum.size())
		{ }

		explicit exact_algorithm(const sbh_data& problem, cxxproperties::Properties& props) :
		exact_algorithm(problem)
		{
			is_position_used = false;
			// is_position_used = props.get<bool>("use-position-info");
		}
	private:
		const size_t sequence_length_;
		const size_t oligo_length_;
		const size_t start_oligo_idx_;
		std::vector<oligo> spectrum_;
		const size_t spectrum_size_;
		std::vector<std::vector<size_t> > graph_;
		bool is_position_used{};
		void generate_graph();
	};
}
