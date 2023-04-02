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
		m_size(size), m_number_of_partcls(4U * N * N * N), m_particles(),
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

	r_point forse_LJ(const r_point& r) inline const noexcept
	{
		return r * (4.0f * m_U_0 / r.abs() / r.abs()) * (12.0f * std::powf(m_sigma / r.abs(), 12U) - 6.0f * std::powf(m_sigma / r.abs(), 6U));
	}


	r_point forse_LJ(const r_point& r1, const r_point& r2) inline const noexcept
	{
		return forse_LJ(r1 - r2);
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

	void init_elementary_cell(const r_point& r_0, float a, float m, float dt);

	void initialize(float m, float dt);

	void initialize_lattice(std::size_t l, float m, float dt);

	void update(float dt);
};