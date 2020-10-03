#include "NelderMead.h"

void NelderMead::checkPoint(std::vector<double>& point)
{
	for (int i = 0; i < this->dimension; i++)
		if (point[i] < lowerBound[i])
			point[i] = lowerBound[i];
		else if (point[i] > upperBound[i])
			point[i] = upperBound[i];
}

NelderMead::NelderMead() {
}

NelderMead::NelderMead(const std::function<double(std::vector<double>)>& costFunction,
	const std::vector<double>& lowerBound, const std::vector<double>& upperBound)
{
	this->costFunction = costFunction;
	this->lowerBound = lowerBound;
	this->upperBound = upperBound;
	this->dimension = static_cast<int>(this->lowerBound.size());
}

NelderMead::~NelderMead() {
}

std::pair<std::vector<double>, double> NelderMead::getResult(bool useTimeLimitCriterion, double stopCriterion)
{	
	std::vector<std::pair<std::vector<double>, double>> simplex;
	for (int i = 0; i <= this->dimension; i++) {
		std::vector<double> tmp;
		for (int j = 0; j < dimension; j++) {
			double low = this->lowerBound[j];
			double high = this->upperBound[j];
			tmp.push_back(low + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (high - low))));
		}
		simplex.push_back({ tmp, 0.0 });
	}

	for (int i = 0; i <= this->dimension; i++)
		simplex[i].second = this->costFunction(simplex[i].first);

	auto comparator = [](const std::pair<std::vector<double>, double>& lhs,
		const std::pair<std::vector<double>, double>& rhs) {
		return lhs.second > rhs.second;
	};

	double f_h, f_g, f_l;
	std::vector<double> x_h, x_g, x_l;
	bool work = true;

	auto tic = std::chrono::high_resolution_clock::now();
	while (work) {

		this->iterations++;
		// проверка условия выхода
		double max = 0.0;
		for (int i = 1; i <= this->dimension; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < this->dimension; j++)
				sum += (simplex[i].first[j] - simplex[0].first[j]) * (simplex[i].first[j] - simplex[0].first[j]);
			sum = std::sqrt(sum);
			if (sum > max)
				max = sum;
		}
		if (!useTimeLimitCriterion && max < stopCriterion)
			work = false;
		if (useTimeLimitCriterion) {
			auto tac = std::chrono::high_resolution_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(tac - tic).count();
			if (elapsed > stopCriterion)
				work = false;
		}

		// сортировка
		std::sort(simplex.begin(), simplex.end(), comparator);

		int idx_h = 0;
		f_h = simplex[idx_h].second;
		x_h = simplex[idx_h].first;

		int idx_g = -1;
		for (int i = 1; i <= dimension; i++)
			if (simplex[i].second < simplex[0].second) {
				idx_g = i;
				break;
			}
		if (idx_g == -1)
			idx_g = 1;
		f_g = simplex[idx_g].second;
		x_g = simplex[idx_g].first;

		int idx_l = this->dimension;
		f_l = simplex[idx_l].second;
		x_l = simplex[idx_l].first;

		// центр тяжести
		std::vector<double> x_c;
		for (int j = 0; j < this->dimension; j++) {
			double sum = 0.0;
			for (int i = 1; i <= this->dimension; i++)
				sum += simplex[i].first[j];
			x_c.push_back(sum / this->dimension);
		}

		// отражение
		std::vector<double> x_r;
		for (int i = 0; i < this->dimension; i++)
			x_r.push_back((1 + this->alpha) * x_c[i] - this->alpha * x_h[i]);
		this->checkPoint(x_r);
		double f_r = this->costFunction(x_r);

		if (f_r < f_l) {
			std::vector<double> x_e;
			for (int i = 0; i < this->dimension; i++)
				x_e.push_back((1 - this->gamma) * x_c[i] + this->gamma * x_r[i]);
			this->checkPoint(x_e);
			double f_e = this->costFunction(x_e);

			if (f_e < f_l) {
				simplex[idx_h] = { x_e, f_e };
				continue;
			}
			else {
				simplex[idx_h] = { x_r, f_r };
				continue;
			}
		}

		if (f_l < f_r && f_r < f_g) {
			simplex[idx_h] = { x_r, f_r };
			continue;
		}

		if (f_h > f_r && f_r > f_g) {
			auto tmp = std::make_pair(x_r, f_r);
			x_r = x_h;
			f_r = f_h;
			x_h = tmp.first;
			f_h = tmp.second;
		}

		// сжатие
		std::vector<double> x_s;
		for (int i = 0; i < this->dimension; i++)
			x_s.push_back(this->beta * x_h[i] + (1 - this->beta) * x_c[i]);
		this->checkPoint(x_s);
		double f_s = this->costFunction(x_s);

		if (f_s < f_h) {
			simplex[idx_h] = { x_s, f_s };
			continue;
		}
		else {
			for (int j = 0; j < this->dimension; j++) // сжимаем всё, кроме последней
				for (int i = 0; i < this->dimension; i++)
					simplex[j].first[i] = x_l[i] + (simplex[j].first[i] - x_l[i]) / 2.;
		}
	}
	return { x_h, f_h };
}