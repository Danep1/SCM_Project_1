#include "Processor.hpp"

void Processor::start()
{
	std::ofstream file_dots;
	std::ofstream file_energy;

	file_dots.open("data.dump");
	file_energy.open("energy.txt");

	//write_current_system(file_dots, 0.0);
	//write_energy(file_energy, 0.0);


	for (auto step = 0U; step < m_N_steps / m_N_update; ++step)
	{
		for (auto update_step = 0U; update_step < m_N_update; ++update_step)
		{
			m_cell.update(dt);
		}

		write_current_system(file_dots, step * m_N_update * dt);
		write_energy(file_energy, step * m_N_update * dt);
	}

	file_dots.close();
	file_energy.close();
}

void Processor::write_energy(std::ofstream& fstream, number_t t)
{
	fstream << t << "\t" << m_cell.get_T() << '\t' << m_cell.get_U() << '\t' << m_cell.get_E() << std::endl;
}

void Processor::write_current_system(std::ofstream& fstream, number_t t)
{
	auto partcls = m_cell.get_particles_ptr();
	fstream << "ITEM: TIMESTEP\n" << t << "\nITEM: NUMBER OF ATOMS\n" << std::size(*partcls) << "\nITEM: BOX BOUNDS pp pp pp\n" << 0.0 << " " << m_length << "\n" << 0.0 << " " << m_width << "\n" << 0.0 << " " << m_height << "\nITEM: ATOMS id x y z E\n";
	auto i = 0U;
	for (auto prtc = std::begin(*partcls); prtc != std::end(*partcls); ++prtc)
	{
		fstream << i << " " << (*prtc).get()->get_pos() << " " << m_cell.get_E() << "\n";
		++i;
	}
}