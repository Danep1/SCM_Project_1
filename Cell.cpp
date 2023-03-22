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

void Cell:: initialize_lattice(std::size_t l, float m, float dt)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution < float > maxwell_dist{ 0.0f, m_v_max };
	std::uniform_real_distribution < float > un_dist;
	r_point r_init;
	r_point v_init;
	auto a = m_size.x() / l;
	for (auto i = 0U; i < l; ++i)
	{
		for (auto j = 0U; j < l; ++j)
		{
			for (auto k = 0U; k < l; ++k)
			{
				r_init = r_point(a * (k + 0.5f), a * (j + 0.5f), a * (i + 0.5f));
				v_init = r_point(maxwell_dist(gen), maxwell_dist(gen), maxwell_dist(gen));
				m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, r_init, r_init + v_init * dt, v_init }));
			}
		}
	}
}

void Cell::update(float dt)
{
	//auto r_min = 2.0f * R_skin;
	m_T = 0.0f;
	m_U = 0.0f;

	r_point a;
	for (auto i = 0U; i < m_number_of_partcls; ++i)
	{
		auto cur_part_ptr = m_particles[i].get();
		cur_part_ptr->save_cur_pos_as_prev();

		a = r_point(0.0f, 0.0f, 0.0f);

		for (auto j = 0U; j < m_number_of_partcls; ++j)
		{
			if (distance(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos()) < m_R_cut && i != j)
			{
				m_U += potential_LJ(m_particles[i], m_particles[j]) / 2.0f;
				a = a + forse_LJ(m_particles[i], m_particles[j]);
			}

			/*if (distance(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos()) < R_cut && i != j)
			{
				m_U += potential_garmonic(m_particles[i], m_particles[j]) / 2.0f;
				a = a + forse_garmonic(m_particles[i], m_particles[j]);
			}*/

			for (auto& trans : m_period_cond_trans)
			{
				if (distance(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos() + trans) < m_R_cut)
				{
					m_U += potential_LJ(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos() + trans) / 2.0f;
					a = a + forse_LJ(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos() + trans) * (1 / cur_part_ptr->get_mass());
				}
				if (distance(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos() - trans) < m_R_cut)
				{
					m_U += potential_LJ(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos() - trans) / 2.0f;
					a = a + forse_LJ(cur_part_ptr->get_pos(), m_particles[j].get()->get_pos() - trans) * (1 / cur_part_ptr->get_mass());
				}
			}
		}
		m_T += cur_part_ptr->get_T();

		cur_part_ptr->accelorate_with(a * dt);
		cur_part_ptr->move_with(cur_part_ptr->get_v() * dt + a * dt * (dt / 2));

		m_E = m_T + m_U;

		if (cur_part_ptr->get_pos().x() < 0.0f)
		{
			cur_part_ptr->move_with(r_point(m_size.x(), 0.0f, 0.0f));
		}
		if (cur_part_ptr->get_pos().x() > m_size.x())
		{
			cur_part_ptr->move_with(r_point(-m_size.x(), 0.0f, 0.0f));
		}
		if (cur_part_ptr->get_pos().y() < 0.0f)
		{
			cur_part_ptr->move_with(r_point(0.0f, m_size.y(), 0.0f));
		}
		if (cur_part_ptr->get_pos().y() > m_size.y())
		{
			cur_part_ptr->move_with(r_point(0.0f, -m_size.y(), 0.0f));
		}
		if (cur_part_ptr->get_pos().z() < 0.0f)
		{
			cur_part_ptr->move_with(r_point(0.0f, 0.0f, m_size.z()));
		}
		if (cur_part_ptr->get_pos().z() > m_size.z())
		{
			cur_part_ptr->move_with(r_point(0.0f, 0.0f, -m_size.z()));
		}
	}
}