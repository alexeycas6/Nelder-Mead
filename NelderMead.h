#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#include <vector>
#include <functional>
#include <algorithm>
#include <chrono>
#include <cmath>

/**
 * \brief Метод оптимизации Нелдер-Мида
 *
 * Эффективен для размерности < 6
 */
class NelderMead
{
private:

	/** Целевая функция */
	std::function<double(std::vector<double>)> costFunction;

	/** Нижняя граница параметров */
	std::vector<double> lowerBound;

	/** Верхняя граница параметров*/
	std::vector<double> upperBound;	

	/** Размерность пространства */
	int dimension;	

	/** Количество итераций */
	int iterations;

	/** Проверка точки на нахождение в пределах границ */
	void checkPoint(std::vector<double>& point);

public:

	/** Коэффициент отражения */
	double alpha{ 1.0 };

	/** Коэффициент сжатия */
	double beta{ 0.5 };

	/** Коэффициент растяжения */
	double gamma{ 2.0 };

	/** Конструктор по умолчанию */
	NelderMead();

	/**
	 * Конструктор
	 * @param costFunction целевая функция
	 * @param lowerBound нижняя граница параметров
	 * @param upperBound верхняя граница параметров
	 */
	NelderMead(const std::function<double(std::vector<double>)>& costFunction,
		const std::vector<double>& lowerBound, const std::vector<double>& upperBound);

	/** Деструктор */
	~NelderMead();

	/**
	 * Расчётный процесс
	 * @param useTimeLimitCriterion если True, то используется ограничение времени работы функции в секундах, иначе по близости точек симплекса
	 * @param stopCriterion критерий остановки (время работы [с] / близость точек симплекса)
	 * @return пару из точки минимума и значения в этой точке
	 */
	std::pair<std::vector<double>, double> getResult(bool useTimeLimitCriterion, double stopCriterion);
};

#endif

