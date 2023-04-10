#pragma once
#include "Cell.hpp"
#include "Functions.hpp"
#include <fstream>
#include <iostream>

class Processor
{
private:
	const float m_length = 1.0f;
	const float m_width = 1.0f;
	const float m_height = 1.0f;
	
	const float m_v_max = 0.01f;

	const std::size_t m_N_particls_in_row = 5U;

	const float m_sigma = m_length / m_N_particls_in_row / 1.095f / 2.0f;
	const float m_R_cut = m_sigma * 2.5f;
	const float m_U_0 = 1.0f;

	const float t = 1.0e+0f;
	const float dt = 1.0e-2f;
	const std::size_t m_N_steps = t / dt;
	const std::size_t m_N_update = 1U;

	Cell m_cell;

public:
	Processor() : m_cell(m_sigma, m_R_cut, m_U_0, m_v_max, r_point(m_length, m_width, m_height), m_N_particls_in_row, dt) {}

	~Processor() noexcept = default;

	void start();

	void write_energy(std::ofstream& fstream, float t);

	void write_current_system(std::ofstream& fstream, float t);
};