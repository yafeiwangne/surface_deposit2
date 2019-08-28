// operate on C_R or C_O, in coupling term they have the same form.
#include "FluxBCudot.h"
#include "Function.h"

registerMooseObject("SurfaceDeposit2App", FluxBCudot);

template <>
InputParameters
validParams<FluxBCudot>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addCoupledVar("coupled_var", "	Choose the variable you want to couple");
  params.addParam<Real>("K_O", 1.0, "standard reaction rate");
  params.addParam<Real>("sigma", 0.01, "deposition layer");
  params.addRequiredParam<FunctionName>("Exp_1", "Exp_1 ");
  params.addRequiredParam<FunctionName>("Exp_2", "Exp_2 ");

  return params;
}

FluxBCudot::FluxBCudot(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _couple_var(coupledValue("coupled_var")),
    _couple_var_offjac(coupled("coupled_var")),
    _u_dot(dot()),
    _du_dot_du(dotDu()),
    _K_O(getParam<Real>("K_O")),
    _sigma(getParam<Real>("sigma")),
    _func_1(getFunction("Exp_1")),
    _func_2(getFunction("Exp_2"))
{
}

Real
FluxBCudot::computeQpResidual()
{
    return  _test[_i][_qp] * (_u_dot[_qp]*_sigma + _K_O*(_u[_qp]*_func_1.value(_t, _q_point[_qp])-_couple_var[_qp] * _func_2.value(_t, _q_point[_qp])));
}

Real
FluxBCudot::computeQpJacobian()
{
  return _test[_i][_qp] * (_du_dot_du[_qp]+ _K_O*_func_1.value(_t, _q_point[_qp])) * _phi[_j][_qp];
}

Real
FluxBCudot::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _couple_var_offjac)
    return _test[_i][_qp] * (-_K_O*_func_2.value(_t, _q_point[_qp])*_phi[_j][_qp]);
  else
    return 0;
}



