#ifndef FLUXBCUDOT_H
#define FLUXBCUDOT_H

#include "IntegratedBC.h"

// Forward Declarations
class FluxBCudot;
class Function;

template <>
InputParameters validParams<FluxBCudot>();

/**
 * Implements Neumann BC where grad(u)=udot-something on the boundary.
 */

class FluxBCudot : public IntegratedBC

{
public:

  FluxBCudot(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  const VariableValue & _couple_var;
  const unsigned _couple_var_offjac;
  const VariableValue & _u_dot;
  const VariableValue & _du_dot_du;
  Real _K_O;
  Real _sigma;
  Function & _func_1;
  Function & _func_2;
};

#endif
