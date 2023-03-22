#include "Cell.hpp"

void Cell::initialize(float m, float dt)
{
	if (m_number_of_partcls == 2U)
	{
		auto v_init = r_point(0.0f, 0.0f, 0.0f);
		auto c = m_size * 0.5f;
		m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, c + r_point(9.0f, 0.0f, 0.0f), c, v_init}));
		m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, c - r_point(9.0f, 0.0f, 0.0f), c, v_init }));
	}
	else
	{ 
		std::random_device rd;
		std::mt19937 gen(rd()); 
		std::normal_distribution < float > maxwell_dist{ 0.0f, m_v_max };
		std::uniform_real_distribution < float > un_dist;
		for (auto i = 0U; i < m_number_of_partcls; ++i)
		{
			auto r_init = r_point(un_dist(gen) * m_size.x(), un_dist(gen) * m_size.y(), un_dist(gen) * m_size.z());
			auto v_init = r_point(maxwell_dist(gen), maxwell_dist(gen), maxwell_dist(gen));
			m_particles.emplace_back(std::make_shared<Particle> (Particle{ m, 0.0f, r_init, r_init + v_init * dt, v_init }));
		}
	}
}