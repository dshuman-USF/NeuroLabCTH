#include <vector>

struct GLMFit
{
  std::vector<double> mu;
  std::vector<double> coef;
};

GLMFit
glmfitpx (std::vector<double>& x, std::vector<double>& y,
         std::vector<double>& offset);
