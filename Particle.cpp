#include "Particle.hpp"

float distance(const r_point& v1, const r_point& v2)
{
	return  std::sqrtf(std::powf(v1.x() - v2.x(), 2U) + std::powf(v1.y() - v2.y(), 2U) + std::powf(v1.z() - v2.z(), 2U));
}

void Particle::move_with(const r_point& dr) noexcept
{
	m_pos = m_pos + dr;
}

void Particle::accelorate_with(const r_point& dv) noexcept
{
	m_v = m_v + dv;
}