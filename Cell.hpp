#pragma once
#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include "Functions.hpp"

class Cell
{
public:
	using particle_t = std::shared_ptr < Particle >;
	using array_t = std::vector < particle_t >;

private:
	const float m_sigma;
	const float m_R_cut;
	const float m_U_0;
	const float m_U_cut; // = 4.0f * U_0 * (std::powf(sigma / R_cut, 12U) - std::powf(sigma / R_cut, 6U));
	const float m_v_max;

	const std::size_t m_number_of_partcls;
	const r_point m_size;

	float m_E = 0.0f;
	float m_T = 0.0f;
	float m_U = 0.0f;

	const std::vector < r_point > m_period_cond_trans;
	array_t m_particles;

public:
	explicit Cell(float sigma, float R_cut, float U_0, float v_max, const r_point& size, std::size_t N, float dt) noexcept:
		m_sigma(sigma), m_R_cut(R_cut), m_U_0(U_0), m_U_cut(4.0f * U_0 * (std::powf(sigma / R_cut, 12U) - std::powf(sigma / R_cut, 6U))),
		m_v_max(v_max),
		m_size(size), m_number_of_partcls(N * N * N), m_particles(),
		m_period_cond_trans({r_point(size.x(), 0.0f, 0.0f), r_point(0.0f, size.y(), 0.0f), r_point(0.0f, 0.0f,  size.z()), r_point(size.x(), size.y(), 0.0f), r_point(size.x(), -size.y(), 0.0f), 
			r_point(size.x(), 0.0f, size.z()),  r_point(size.x(), 0.0f, -size.z()), r_point(0.0f, size.y(), size.z()), r_point(0.0f, size.y(), -size.z()), 
			r_point(size.x(), size.y(), size.z()), r_point(size.x(), size.y(), -size.z()), r_point(size.x(), -size.y(), size.z()), r_point(-size.x(), size.y(), size.z())})
	{
		initialize_lattice(N, 1.0f, dt);
	}

	~Cell() noexcept = default;

public:
	auto get_particles_ptr() inline const noexcept
	{
		return std::make_shared <array_t> (m_particles);
	}

	auto  get_E() inline const noexcept
	{
		return m_E;
	}

	auto  get_T() inline const noexcept
	{
		return m_T;
	}

	auto  get_U() inline const noexcept
	{
		return m_U;
	}

	float potential_LJ(const r_point& r) inline const noexcept
	{
		return 4.0f * m_U_0 * (std::powf(m_sigma / r.abs(), 12U) - std::powf(m_sigma / r.abs(), 6U)) - m_U_cut;
	}

	float potential_LJ(const r_point& r1, const r_point& r2) inline const noexcept
	{
		return potential_LJ(r2 - r1);
	}

	float potential_LJ(const particle_t& p1, const particle_t& p2) inline const noexcept
	{
		return potential_LJ(p2.get()->get_pos() - p1.get()->get_pos());
	}

	r_point forse_LJ(const r_point& r1, const r_point& r2) inline const noexcept
	{
		auto r = r1 - r2;
		return r * (4.0f * m_U_0 / r.abs() / r.abs()) * (12.0f * std::powf(m_sigma / r.abs(), 12U) - 6.0f * std::powf(m_sigma / r.abs(), 6U));
	}

	r_point forse_LJ(const particle_t& p1, const particle_t& p2) inline const noexcept
	{
		return forse_LJ(p1.get()->get_pos(), p2.get()->get_pos());
	}

	//float potential_garmonic(const r_point& r) inline const noexcept
	//{
	//	return std::powf(r.abs() - 0.1f, 2U) / 2;
	//}

	//float potential_garmonic(const r_point& r1, const r_point& r2) inline const noexcept
	//{
	//	return potential_garmonic(r2 - r1);
	//}

	//float potential_garmonic(const particle_t& p1, const particle_t& p2) inline const noexcept
	//{
	//	return potential_garmonic(p2.get()->get_pos() - p1.get()->get_pos());
	//}

	//r_point forse_garmonic(const r_point& r1, const r_point& r2) inline const noexcept
	//{
	//	auto r = r2 - r1;
	//	return r * ((r.abs() - 0.1f) / r.abs());
	//}

	//r_point forse_garmonic(const particle_t& p1, const particle_t& p2) inline const noexcept
	//{
	//	return forse_garmonic(p1.get()->get_pos(), p2.get()->get_pos());
	//}


	void initialize(float m, float dt);

	void initialize_lattice(std::size_t l, float m, float dt)
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
					r_init = r_point( a * (k + 0.5f), a * (j + 0.5f), a * (i + 0.5f));
					v_init = r_point(maxwell_dist(gen), maxwell_dist(gen), maxwell_dist(gen));
					m_particles.emplace_back(std::make_shared<Particle>(Particle{ m, 0.0f, r_init, r_init + v_init * dt, v_init }));
				}
			}
		}
	}

	void update(float dt)
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
				cur_part_ptr->move_with( r_point(m_size.x(), 0.0f, 0.0f) );
			}
			if (cur_part_ptr->get_pos().x() > m_size.x())
			{
				cur_part_ptr->move_with( r_point(-m_size.x(), 0.0f, 0.0f) );
			}
			if (cur_part_ptr->get_pos().y() < 0.0f)
			{
				cur_part_ptr->move_with(r_point(0.0f, m_size.y(), 0.0f));
			}
			if (cur_part_ptr->get_pos().y() > m_size.y())
			{
				cur_part_ptr->move_with(r_point(0.0f, - m_size.y(), 0.0f));
			}
			if (cur_part_ptr->get_pos().z() < 0.0f)
			{
				cur_part_ptr->move_with(r_point(0.0f, 0.0f, m_size.z()));
			}
			if (cur_part_ptr->get_pos().z() > m_size.z())
			{
				cur_part_ptr->move_with(r_point(0.0f, 0.0f, - m_size.z()));
			}
		}
	}
};