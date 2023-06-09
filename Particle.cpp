#include "Particle.hpp"

number_t distance(const r_point& v1, const r_point& v2)
{
	return  std::sqrtf(std::powf(v1.x() - v2.x(), 2U) + std::powf(v1.y() - v2.y(), 2U) + std::powf(v1.z() - v2.z(), 2U));
}

void Particle::update(number_t dt) noexcept
{
	this->save_cur_pos_as_prev();
	m_v = m_v + m_a * dt; // ������ � ����� �������, ����� ����� ���� �������
	m_pos = m_pos + m_v * dt + m_a * dt * (dt / 2.0);		
	m_a = r_point(0.0, 0.0, 0.0);
}


void Particle::move_with(const r_point& dr) noexcept
{
	m_pos = m_pos + dr;
}