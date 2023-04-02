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
		r_init = r_0 + (*pos) * (a / 2.0f);
		v_init = r_point(v_dist(gen), v_dist(gen), v_dist(gen)) * m_v_max;
		m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, r_init, r_init + v_init * dt, v_init }));
	}
}

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
	//auto r_min = 2.0f * R_skin;
	m_T = 0.0f;
	m_U = 0.0f;

	r_point a;
	for (auto i = std::begin(m_particles); i != std::end(m_particles); ++i)
	{
		i->get()->save_cur_pos_as_prev();

		a = r_point(0.0f, 0.0f, 0.0f);
		std::vector <r_point> periodic_positions;
		auto dx = 0.0f;
		auto dy = 0.0f;
		auto dz = 0.0f;

		for (auto j = std::begin(m_particles); j != std::end(m_particles); ++j)
		{
			dx = i->get()->get_pos().x() - j->get()->get_pos().x();
			dy = i->get()->get_pos().y() - j->get()->get_pos().y();
			dz = i->get()->get_pos().z() - j->get()->get_pos().z();
			dx = std::fmin(dx, std::fmin(dx + m_size.x(), dx - m_size.x()));
			dy = std::fmin(dy, std::fmin(dy + m_size.y(), dy - m_size.y()));
			dz = std::fmin(dz, std::fmin(dz + m_size.z(), dz - m_size.z()));
			auto r = r_point(dx, dy, dz);
			if (r.abs() < m_R_cut)
			{
				m_U += potential_LJ(r) / 2.0f;
				a = a + forse_LJ(r) * (1 / i->get()->get_mass());
			}
		}

		m_T += i->get()->get_T();

		i->get()->move_with(i->get()->get_v() * dt + a * dt * (dt / 2));	// r(t + dt) = r(t) + v(t) * dt + a(t) /2 * dt^2
		i->get()->accelorate_with(a * dt);									// v(t + dt) = v(t) + a(t) * dt

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
	}
}