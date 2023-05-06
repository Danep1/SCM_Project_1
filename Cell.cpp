#include "Cell.hpp"

void Cell::init_elementary_cell(const r_point& r_0, float a, float m, float dt)
{
	auto positions = { r_point(0.0f, 0.0f, 0.0f), r_point(0.0f, 1.0f, 1.0f), r_point(1.0f, 0.0f, 1.0f), r_point(1.0f, 1.0f, 0.0f) };
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution < float > v_dist{ 0.0f, 1.0f };
	r_point r_init;
	r_point v_init;
	for (auto pos = std::begin(positions); pos != std::end(positions); ++pos)
	{
		r_init = r_0 + (*pos) * (a * 0.5f);
		v_init = r_point(v_dist(gen), v_dist(gen), v_dist(gen)) * m_v_max;
		m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, r_init, r_init + v_init * dt, v_init }));
		++m_number_of_partcls;
	}
}

void Cell::initialize_dipole(float m, float dt)
{
	auto v_init = r_point(0.0f, 0.0f, 0.0f);
	auto c = m_size * 0.5f;
	m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, c + r_point(0.15f, 0.0f, 0.0f), c, v_init}));
	m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, c - r_point(0.15f, 0.0f, 0.0f), c, v_init }));
	++m_number_of_partcls;
	++m_number_of_partcls;
}

void Cell:: initialize_lattice(std::size_t l, float m, float dt)
{
	auto a = m_size.x() / l;
	auto r_0 = r_point(1.0f, 1.0f, 1.0f) * 0.05f;
	for (auto i = 0U; i < l; ++i)
	{
		for (auto j = 0U; j < l; ++j)
		{
			for (auto k = 0U; k < l; ++k)
			{
				init_elementary_cell(r_0 + r_point(i, j, k) * a, a, m, dt);
			}
		}
	}
}



void Cell::update(float dt)
{
	m_T = 0.0f;
	m_U = 0.0f;
	r_point a;
	for (auto i = std::begin(m_particles); i != std::end(m_particles); ++i)
	{
		i->get()->save_cur_pos_as_prev();
		std::vector <r_point> periodic_positions;

		for (auto j = std::begin(m_particles); j != i; ++j)
		for (auto j = std::begin(m_particles); j != std::end(m_particles); ++j)
		{
		/*	s = { std::fabsf(i->get()->get_pos().x() - j->get()->get_pos().x()), std::fabsf(i->get()->get_pos().x() - j->get()->get_pos().x() - m_size.x()), std::fabsf(i->get()->get_pos().x() - j->get()->get_pos().x() + m_size.x()) };
			dx = *std::min_element(std::begin(s), std::end(s));

			s = { std::fabsf(i->get()->get_pos().y() - j->get()->get_pos().y()), std::fabsf(i->get()->get_pos().y() - j->get()->get_pos().y() - m_size.y()), std::fabsf(i->get()->get_pos().y() - j->get()->get_pos().y() + m_size.y()) };
			dy = *std::min_element(std::begin(s), std::end(s));

			s = { std::fabsf(i->get()->get_pos().z() - j->get()->get_pos().z()), std::fabsf(i->get()->get_pos().z() - j->get()->get_pos().z() - m_size.z()), std::fabsf(i->get()->get_pos().z() - j->get()->get_pos().z() + m_size.z()) };
			dz = *std::min_element(std::begin(s), std::end(s));*/
			//r = r_point(dx, dy, dz);

			r = i->get()->get_pos() - j->get()->get_pos();
			if (r * r < m_R_cut * m_R_cut && i != j)
			{
				m_U += potential_garmonic(r) / 2.0f;
				a = a + forse_garmonic(i->get()->get_pos(), j->get()->get_pos());
			auto min_r = std::min_element(std::begin(periodic_positions), std::end(periodic_positions), [](auto r1, auto r2) {return r1.abs() < r2.abs(); });
			
			m_U += potential_LJ(*min_r);
			i->get()->m_a = i->get()->m_a + forse_LJ(*min_r);
			j->get()->m_a = j->get()->m_a - forse_LJ(*min_r);
			a = a + forse_LJ(*min_r) * (1 / i->get()->get_mass());
		}
	}	
	for (auto i = std::begin(m_particles); i != std::end(m_particles); ++i)
	{
		i->get()->update(dt);
		m_E = m_T + m_U;

		if (i->get()->get_pos().x() < 0.0f)
		{
			i->get()->move_with(r_point(m_size.x(), 0.0f, 0.0f));
		}
		else if (i->get()->get_pos().x() > m_size.x())
		{
			i->get()->move_with(r_point(-m_size.x(), 0.0f, 0.0f));
		}
		if (i->get()->get_pos().y() < 0.0f)
		{
			i->get()->move_with(r_point(0.0f, m_size.y(), 0.0f));
		}
		else if (i->get()->get_pos().y() > m_size.y())
		{
			i->get()->move_with(r_point(0.0f, -m_size.y(), 0.0f));
		}
		else if (i->get()->get_pos().z() < 0.0f)
		{
			i->get()->move_with(r_point(0.0f, 0.0f, m_size.z()));
		}
		else if (i->get()->get_pos().z() > m_size.z())
		{
			i->get()->move_with(r_point(0.0f, 0.0f, -m_size.z()));
		}

		i->get()->m_a = r_point(0.0f, 0.0f, 0.0f);

		m_T += i->get()->get_T();
		m_E = m_T + m_U;
	}
}