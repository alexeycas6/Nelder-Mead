# Nelder-Mead
Realization of Nelder-Mead optimization algorithm
# Usage
The process can end according to two criteria:
+ closeness simlex points
  ```c++
  NelderMead optimization(function, { -4.0, -4.0, -4.0 }, { 4.0, 4.0, 4.0 });
  auto res_by_closeness = optimization.getResult(false, 1.E-5);
  ```
+ work time (in seconds)
  ```c++
  NelderMead optimization(function, { -4.0, -4.0, -4.0 }, { 4.0, 4.0, 4.0 });
  auto res_by_time = optimization.getResult(true, 3.0);
  ```
