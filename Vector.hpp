#pragma once
#include <cmath>
#include <iostream>

using number_t = double;

class r_point
{
private:
	number_t m_x;
	number_t m_y;
	number_t m_z;

public:
	r_point() noexcept = default;

	explicit r_point(number_t x, number_t y, number_t z) noexcept :
		m_x(x), m_y(y), m_z(z)
	{}

	r_point(const r_point &) = default;

	~r_point() noexcept = default;

public:
	const auto x() inline const noexcept
	{
		return m_x;
	}

	const auto y() inline const noexcept
	{
		return m_y;
	}

	const auto z() inline const noexcept	
	{
		return m_z;
	}

	const auto abs() inline const noexcept
	{
		return std::sqrtf((*this) * (* this));
	}

public:
	r_point operator+ (const r_point& v2) const;

	r_point operator- (const r_point& v2) const;

	r_point operator* (number_t a) const;

	number_t operator* (const r_point& v2) const;

	r_point& operator= (const r_point& other) = default;

	friend std::ostream& operator<<(std::ostream& os, const r_point& dt);
};