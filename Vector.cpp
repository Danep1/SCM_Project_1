#include "vector.hpp"

r_point r_point::operator+ (const r_point& v2) const
{
	return r_point(m_x + v2.x(), m_y + v2.y(), m_z + v2.z());
}

r_point r_point::operator- (const r_point& v2) const
{
	return r_point(m_x - v2.x(), m_y - v2.y(), m_z - v2.z());
}

r_point r_point::operator* (number_t a) const
{
	return r_point(m_x * a, m_y * a, m_z * a);
}

number_t r_point::operator* (const r_point& v2) const
{
	return m_x * v2.x() + m_y * v2.y() + m_z * v2.z();
}


std::ostream& operator<< (std::ostream& os, const r_point& v)
{
	os << v.x() << " " << v.y() << " " << v.z();
	return os;
}