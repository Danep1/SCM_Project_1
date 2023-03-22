#pragma once
#include "Cell.hpp"
#include "Functions.hpp"
#include <fstream>
#include <iostream>

class Processor
{
private:
	const float m_length = 100.0f;
	const float m_width = 100.0f;
	const float m_height = 100.0f;
	
	const float m_v_max = 5.0f;

	const std::size_t m_N_particls_in_row = 4U;

	const float m_sigma = m_length / m_N_particls_in_row / 1.1224620f;
	const float m_R_cut = m_sigma * 2.5f;
	const float m_U_0 = 1.0f;

	const float t = 1.0f;
	const float dt = 0.0001f;
	const std::size_t m_N_steps = t / dt;
	const std::size_t m_N_update = 100U;

	Cell m_cell;

public:
	Processor() : m_cell(m_sigma, m_R_cut, m_U_0, m_v_max, r_point(m_length, m_width, m_height), m_N_particls_in_row, dt) {}

	~Processor() noexcept = default;

	void start()
	{
		std::ofstream file_dots;
		std::ofstream file_energy;

		file_dots.open("data.dump");
		file_energy.open("energy.txt");

		write_current_system(file_dots, 0.0f);
		write_energy(file_energy, 0.0f);

		for (auto step = 0U; step < m_N_steps / m_N_update; ++step)
		{
			for (auto update_step = 0U; update_step < m_N_update; ++update_step) 
			{
				m_cell.update(dt);
			}

			write_current_system(file_dots, step *  dt);
			write_energy(file_energy, step *  dt);
		}
		
		file_dots.close();
		file_energy.close();
	}

	void write_energy(std::ofstream& fstream, float t)
	{
		fstream << t << "\t" << m_cell.get_T() << '\t' << m_cell.get_U() << '\t' << m_cell.get_E() << std::endl;
	}

	void write_current_system(std::ofstream& fstream, float t)
	{
		auto partcls = m_cell.get_particles_ptr();
		fstream << "ITEM: TIMESTEP\n" << t << "\nITEM: NUMBER OF ATOMS\n" << m_N_particls_in_row * m_N_particls_in_row * m_N_particls_in_row << "\nITEM: BOX BOUNDS pp pp pp\n" << 0.0f << " " << m_length << "\n" << 0.0f << " " << m_width << "\n" << 0.0f << " " << m_height << "\nITEM: ATOMS id x y z E\n";
		for (auto i = 0U; i < m_N_particls_in_row * m_N_particls_in_row * m_N_particls_in_row; ++i)
		{
			fstream << i << " " << partcls->at(i).get()->get_pos() << " " << m_cell.get_E() << "\n";
		}
	}
};