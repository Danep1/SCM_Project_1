#pragma once
#include "Cell.hpp"
#include "Functions.hpp"
#include <fstream>
#include <iostream>

class Processor
{
private:
	const number_t m_length = 1.0;
	const number_t m_width = 1.0;
	const number_t m_height = 1.0;
	
	const number_t m_v_max = 1.0e+0;

	const std::size_t m_N_particls_in_row = 5U;

	const number_t m_sigma = m_length / m_N_particls_in_row / 1.095 / std::sqrt(2.0);
	const number_t m_R_cut = m_sigma * 2.5;
	const number_t m_U_0 = 1.0;

	const number_t t = 1.0e-1;
	const number_t dt = 1.0e-5;
	const std::size_t m_N_steps = t / dt;
	const std::size_t m_N_update = 100U;

	Cell m_cell;

public:
	Processor() : m_cell(m_sigma, m_R_cut, m_U_0, m_v_max, r_point(m_length, m_width, m_height), m_N_particls_in_row, dt) {}

	~Processor() noexcept = default;

	void start();

	void write_energy(std::ofstream& fstream, number_t t);

	void write_current_system(std::ofstream& fstream, number_t t);
};