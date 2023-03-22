#include "Processor.hpp"
#include "timer.hpp"

int main()
{
	Timer t1("Calculation", true);
	{
		auto pr = Processor();
		pr.start();
	}
}