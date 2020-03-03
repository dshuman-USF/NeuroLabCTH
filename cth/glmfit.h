#include <vector>

struct GLMFit
{
  std::vector<double> mu;
  std::vector<double> coef;
};

GLMFit
glmfitp (std::vector<double>& x, std::vector<double>& y,
         std::vector<double>& offset);
