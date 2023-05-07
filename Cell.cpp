#include "Cell.hpp"

void Cell::init_elementary_cell(const r_point& r_0, number_t a, number_t m, number_t dt)
{
	auto positions = { r_point(0.0, 0.0, 0.0), r_point(0.0, 1.0, 1.0), r_point(1.0, 0.0, 1.0), r_point(1.0, 1.0, 0.0) };
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution < number_t > v_dist{ 0.0, 1.0 };
	r_point r_init;
	r_point v_init;
	for (auto pos = std::begin(positions); pos != std::end(positions); ++pos)
	{
		r_init = r_0 + (*pos) * (a * 0.5);
		v_init = r_point(v_dist(gen), v_dist(gen), v_dist(gen)) * m_v_max;
		m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0, r_init, r_init + v_init * dt, v_init }));
		++m_number_of_partcls;
	}
}

void Cell::initialize_dipole(number_t m, number_t dt)
{
	auto v_init = r_point(0.0, 0.0, 0.0);
	auto c = m_size * 0.5;
	m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0, c + r_point(1.0, 0.0, 0.0) * 0.43, c, v_init }));
	m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0, c - r_point(1.0, 0.0, 0.0) * 0.43, c, v_init }));
	++m_number_of_partcls;
	++m_number_of_partcls;
}

void Cell::initialize_lattice(std::size_t l, number_t m, number_t dt)
{
	auto a = m_size.x() / l;
	auto r_0 = r_point(1.0, 1.0, 1.0) * 0.05;
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


void Cell::update(number_t dt)
{
	m_T = 0.0;
	m_U = 0.0;

	number_t dx, dy, dz;
	std::vector <number_t> s;
	r_point min_r;

	for (auto i = std::begin(m_particles); i != std::end(m_particles); ++i)
	{
		for (auto j = std::begin(m_particles); j != i; ++j)
		{
			dx = i->get()->get_pos().x() - j->get()->get_pos().x();
			s = {dx, dx + m_size.x(), dx - m_size.x()};
			dx = *std::min_element(std::begin(s), std::end(s), [](auto l, auto r) {return std::abs(l) < std::abs(r); });

			dy = i->get()->get_pos().y() - j->get()->get_pos().y();
			s = { dy, dy + m_size.y(), dy - m_size.y() };
			dy = *std::min_element(std::begin(s), std::end(s), [](auto l, auto r) {return std::abs(l) < std::abs(r); });

			dz = i->get()->get_pos().z() - j->get()->get_pos().z();
			s = { dz, dz + m_size.z(), dz - m_size.z() };
			dz = *std::min_element(std::begin(s), std::end(s), [](auto l, auto r) {return std::abs(l) < std::abs(r); });

			min_r = r_point(dx, dy, dz);

			if (min_r * min_r < m_R_cut * m_R_cut)
			{
				m_U += potential_LJ(min_r);
				i->get()->m_a = i->get()->m_a + forse_LJ(min_r);
				j->get()->m_a = j->get()->m_a - forse_LJ(min_r);
			}
		}
	}
	for (auto i = std::begin(m_particles); i != std::end(m_particles); ++i)
	{

		i->get()->update(dt);

		if (i->get()->get_pos().x() < 0.0)
		{
			i->get()->move_with(r_point(m_size.x(), 0.0, 0.0));
		}
		else if (i->get()->get_pos().x() > m_size.x())
		{
			i->get()->move_with(r_point(-m_size.x(), 0.0, 0.0));
		}
		if (i->get()->get_pos().y() < 0.0)
		{
			i->get()->move_with(r_point(0.0, m_size.y(), 0.0));
		}
		else if (i->get()->get_pos().y() > m_size.y())
		{
			i->get()->move_with(r_point(0.0, -m_size.y(), 0.0));
		}
		else if (i->get()->get_pos().z() < 0.0)
		{
			i->get()->move_with(r_point(0.0, 0.0, m_size.z()));
		}
		else if (i->get()->get_pos().z() > m_size.z())
		{
			i->get()->move_with(r_point(0.0, 0.0, -m_size.z()));
		}

		m_T += i->get()->get_T();
		m_E = m_T + m_U;
	}
}