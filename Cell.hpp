#pragma once
#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include "Particle.hpp"

class Cell
{
public:
	using particle_t = std::shared_ptr < Particle >;
	using array_t = std::vector < particle_t >;

private:
	const number_t m_sigma;
	const number_t m_R_cut;
	const number_t m_U_0;
	const number_t m_U_cut = 4.0 * m_U_0 * (std::pow(m_sigma / m_R_cut, 12U) - std::pow(m_sigma / m_R_cut, 6U));
	const number_t m_v_max;

	std::size_t m_number_of_partcls;
	const r_point m_size;

	number_t m_E = 0.0;
	number_t m_T = 0.0;
	number_t m_U = 0.0;

	const std::vector < r_point > m_period_cond_trans;
	array_t m_particles;

public:
	explicit Cell(number_t sigma, number_t R_cut, number_t U_0, number_t v_max, const r_point& size, std::size_t N, number_t dt) noexcept :
		m_sigma(sigma), m_R_cut(R_cut), m_U_0(U_0), m_U_cut(4.0 * U_0 * (std::powf(sigma / R_cut, 12U) - std::powf(sigma / R_cut, 6U))),
		m_v_max(v_max),
		m_size(size), m_number_of_partcls(0U), m_particles(),
		m_period_cond_trans({ r_point(size.x(), 0.0, 0.0), r_point(0.0, size.y(), 0.0), r_point(0.0, 0.0,  size.z()), r_point(size.x(), size.y(), 0.0), r_point(size.x(), -size.y(), 0.0),
			r_point(size.x(), 0.0, size.z()),  r_point(size.x(), 0.0, -size.z()), r_point(0.0, size.y(), size.z()), r_point(0.0, size.y(), -size.z()),
			r_point(size.x(), size.y(), size.z()), r_point(size.x(), size.y(), -size.z()), r_point(size.x(), -size.y(), size.z()), r_point(-size.x(), size.y(), size.z()) })
	{
		initialize_lattice(N, 1.0, dt);
	}

	~Cell() noexcept = default;

public:
	auto get_particles_ptr() inline const noexcept
	{
		return std::make_shared <array_t>(m_particles);
	}

	auto get_particles_begin() inline const noexcept
	{
		return std::begin(m_particles);
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

	number_t potential_LJ(const r_point& r) inline const noexcept
	{
		return 4.0 * m_U_0 * (std::powf(m_sigma / r.abs(), 12U) - std::powf(m_sigma / r.abs(), 6U)) - m_U_cut;
	}

	number_t potential_LJ(const r_point& r1, const r_point& r2) inline const noexcept
	{
		return potential_LJ(r2 - r1);
	}

	number_t potential_LJ(const particle_t& p1, const particle_t& p2) inline const noexcept
	{
		return potential_LJ(p2.get()->get_pos() - p1.get()->get_pos());
	}

	r_point forse_LJ(const r_point& r) inline const noexcept
	{
		return r * (4.0 * m_U_0 / r.abs() / r.abs()) * (12.0 * std::pow(m_sigma / r.abs(), 12U) - 6.0 * std::pow(m_sigma / r.abs(), 6U));
	}

	r_point forse_LJ(const r_point& r1, const r_point& r2) inline const noexcept
	{
		return forse_LJ(r1 - r2);
	}

	r_point forse_LJ(const particle_t& p1, const particle_t& p2) inline const noexcept
	{
		return forse_LJ(p1.get()->get_pos(), p2.get()->get_pos());
	}

	number_t potential_garmonic(const r_point& r) inline const noexcept
	{
		return std::powf(r.abs() - 0.25, 2U) / 2;
	}

	number_t potential_garmonic(const r_point& r1, const r_point& r2) inline const noexcept
	{
		return potential_garmonic(r2 - r1);
	}

	number_t potential_garmonic(const particle_t& p1, const particle_t& p2) inline const noexcept
	{
		return potential_garmonic(p2.get()->get_pos() - p1.get()->get_pos());
	}

	r_point forse_garmonic(const r_point& r, number_t r0) inline const noexcept
	{
		return r * ((r.abs() - r0) / r.abs());
	}

	r_point forse_garmonic(const r_point& r1, const r_point& r2, number_t r0) inline const noexcept
	{
		return forse_garmonic(r2 - r1, r0);
	}

	r_point forse_garmonic(const particle_t& p1, const particle_t& p2, number_t r0) inline const noexcept
	{
		return forse_garmonic(p1.get()->get_pos(), p2.get()->get_pos(), r0);
	}

	void init_elementary_cell(const r_point& r_0, number_t a, number_t m, number_t dt);

	void initialize_dipole(number_t m, number_t dt);

	void initialize_lattice(std::size_t l, number_t m, number_t dt);

	void update(number_t dt);
};