#pragma once

#include <utility>
#include <vector>
#include <string>
#include <iomanip>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"


enum nucleotides {
	adenine, cytosine, guanine, thymine
};

const char nucleotides_array[4] = { 'A', 'C', 'G', 'T' };

struct oligo {
	oligo(const size_t index, const int pos_l, const int pos_h, const uint32_t seq)
		: index(index),
		  posL(pos_l),
		  posH(pos_h),
		  seq(seq) {}

	oligo()
		: index(0),
		  posL(0),
		  posH(0),
		  seq(0u) {}
	

	size_t index;
	size_t posL;
	size_t posH;
	uint32_t seq;

	bool operator <(const oligo& rhs) const
	{
		return seq < rhs.seq;
	}

	friend bool operator <(const uint32_t seq, const oligo& rhs) 
	{
		return seq < rhs.seq;
	}
};

struct sbh_data
{
	sbh_data(const size_t oligo_length, const size_t sequence_length, const size_t start_oligo_idx,
	         std::vector<oligo> spectrum)
		: oligo_length(oligo_length),
		  sequence_length(sequence_length),
		  start_oligo_idx(start_oligo_idx),
		  spectrum(std::move(spectrum)) {}

	const size_t oligo_length;
	const size_t sequence_length;
	const size_t start_oligo_idx;
	std::vector <oligo> spectrum;
};

inline uint32_t encode(const std::string& oligo) {
	assert(oligo.size() < 16);
	auto result = 0;
	for (const char nucleotide : oligo) {
		result <<= 2;
		switch (nucleotide) {
		case 'A': result |= adenine; break;
		case 'C': result |= cytosine; break;
		case 'G': result |= guanine; break;
		case 'T': result |= thymine; break;
		default: break;
		}
	}
	return result;
}

inline std::string decode_n_last(const uint32_t oligo, const size_t oligo_length) {
	std::string result;
	result.resize(oligo_length);
	for (size_t i = 0; i < oligo_length; ++i)
		result[oligo_length - i - 1] = nucleotides_array[(oligo >> (2 * i)) & 3];
	return result;
}

inline size_t calc_overlap(const uint32_t lhs, const uint32_t rhs, const size_t oligo_length) {
	const uint32_t mask = (1 << (2 * oligo_length)) - 1; //set 2*k last bits
	auto i = 2;
	auto overlap = oligo_length;
	while (overlap-- && (lhs ^ (rhs >> i)) & mask >> i) i += 2;
	return overlap;
}


inline struct sbh_data load_problem_hackerrank(std::istream& is) {

	std::string dummy, start;
	std::vector<oligo> spectrum;
	int sequence_length, oligo_length;
	is >> dummy >> sequence_length;
	is >> dummy >> start;
	is >> dummy >> oligo_length;
	is >> std::ws;

	size_t idx = 0;
	while (getline(is, dummy)) {
		spectrum.emplace_back(idx, 0, 0, encode(dummy));
		idx++;
	}

	const struct oligo start_oligo(0, 0, 0, encode(start));

	const auto start_it = std::lower_bound(spectrum.begin(), spectrum.end(), start_oligo);
	const size_t start_oligo_idx = std::distance(spectrum.begin(), start_it);

	return sbh_data(oligo_length, sequence_length, start_oligo_idx, spectrum);
}

inline struct sbh_data load_problem_xml(std::istream& is) {
	//parse XML unsigned
	std::vector<char> buffer;

	is.seekg(0, std::ios::end);
	buffer.reserve(is.tellg());
	is.seekg(0, std::ios::beg);
	buffer.assign((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char	>());
	char* pt = buffer.data();

	rapidxml::xml_document<>doc;
	doc.parse<0>(pt);

	//get length of final sequence
	rapidxml::xml_node<>* node = doc.first_node();
	rapidxml::xml_attribute<>* attr = node->first_attribute("length");
	const size_t sequence_length = strtol(attr->value(), nullptr, 10);

	//get start oligo
	attr = attr->next_attribute("start");
	const uint32_t start_oligo_seq = encode(std::string(attr->value()));
	const oligo start_oligo{0, 0,0,start_oligo_seq};

	//get length of each oligo (k)
	node = node->first_node("probe");
	attr = node->first_attribute("pattern");
	const size_t oligo_length = attr->value_size();
	assert(oligo_length < 16);

	const size_t spectrum_size = count_children(node);
	std::vector<oligo> spectrum(spectrum_size);


	size_t idx = 0;
	for (node = node->first_node("cell"); node != nullptr; node = node->next_sibling()) {
		spectrum[idx].index = idx;
		attr = node->first_attribute("posL");
		spectrum[idx].posL = strtol(attr->value(), nullptr, 10);

		attr = attr->next_attribute("posH");
		spectrum[idx].posH = strtol(attr->value(), nullptr, 10);

		spectrum[idx].seq = encode(std::string(node->value()));
		idx++;
	}

	const auto start_it = std::lower_bound(spectrum.begin(), spectrum.end(), start_oligo);
	const size_t start_oligo_idx = std::distance(spectrum.begin(), start_it);
	return sbh_data(oligo_length, sequence_length, start_oligo_idx, spectrum);
}
