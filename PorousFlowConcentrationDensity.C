#include "PorousFlowConcentrationDensity.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("PorousFlowApp", PorousFlowConcentrationDensity);
registerMooseObject("PorousFlowApp", ADPorousFlowConcentrationDensity);

template <bool is_ad>
InputParameters
PorousFlowConcentrationDensityTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowFluidPropertiesBaseTempl<is_ad>::validParams();
  params.addRequiredParam<UserObjectName>("fp", "The name of the user object for fluid properties");
  params.addCoupledVar("concentration", 0, "The mass fraction variable");
  params.addClassDescription("This Material calculates fluid properties at the quadpoints or nodes "
                             "for a single component fluid");
  return params;
}

template <bool is_ad>
PorousFlowConcentrationDensityTempl<is_ad>::PorousFlowConcentrationDensityTempl(
    const InputParameters & parameters)
  : PorousFlowFluidPropertiesBaseTempl<is_ad>(parameters),
    _fp(this->template getUserObject<SinglePhaseFluidProperties>("fp")),
    _is_concentration_nodal(isCoupled("concentration") ? getFieldVar("concentration", 0)->isNodal() : false),
    _concentration(_nodal_material && _is_concentration_nodal ? this->template coupledGenericDofValue<is_ad>("concentration")
                                      : this->template coupledGenericValue<is_ad>("concentration"))
{
}

template <bool is_ad>
void
PorousFlowConcentrationDensityTempl<is_ad>::initQpStatefulProperties()
{
  if (_compute_rho_mu)
  {
    (*_density)[_qp] = 1000 + _concentration[_qp] * 0.0001;

    (*_viscosity)[_qp] = _fp.mu_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                         _temperature[_qp] + _t_c2k) /
                         _pressure_to_Pascals / _time_to_seconds;
  }

  if (_compute_internal_energy)
    (*_internal_energy)[_qp] = _fp.e_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                              _temperature[_qp] + _t_c2k);

  if (_compute_enthalpy)
    (*_enthalpy)[_qp] = _fp.h_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                       _temperature[_qp] + _t_c2k);
}

template <bool is_ad>
void
PorousFlowConcentrationDensityTempl<is_ad>::computeQpProperties()
{
  if (_compute_rho_mu)
  {
    if (is_ad)
    {
      GenericReal<is_ad> rho, mu;
      _fp.rho_mu_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                          _temperature[_qp] + _t_c2k,
                          rho,
                          mu);

      (*_density)[_qp] = 1000 + _concentration[_qp] * 0.0001;
      (*_viscosity)[_qp] = mu / _pressure_to_Pascals / _time_to_seconds;
    }
    else
    {
      // Density and viscosity, and derivatives wrt pressure and temperature
      Real rho, drho_dp, drho_dT, mu, dmu_dp, dmu_dT;
      _fp.rho_mu_from_p_T(MetaPhysicL::raw_value(_porepressure[_qp][_phase_num]) *
                              _pressure_to_Pascals,
                          MetaPhysicL::raw_value(_temperature[_qp]) + _t_c2k,
                          rho,
                          drho_dp,
                          drho_dT,
                          mu,
                          dmu_dp,
                          dmu_dT);
      (*_density)[_qp] = 1000 + _concentration[_qp] * 0.0001;
      (*_ddensity_dp)[_qp] = 0;
      (*_ddensity_dT)[_qp] = drho_dT;
      (*_viscosity)[_qp] = mu / _pressure_to_Pascals / _time_to_seconds;
      (*_dviscosity_dp)[_qp] = dmu_dp / _time_to_seconds;
      (*_dviscosity_dT)[_qp] = dmu_dT / _pressure_to_Pascals / _time_to_seconds;
    }
  }

  if (_compute_internal_energy)
  {
    if (is_ad)
      (*_internal_energy)[_qp] = _fp.e_from_p_T(
          _porepressure[_qp][_phase_num] * _pressure_to_Pascals, _temperature[_qp] + _t_c2k);
    else
    {
      // Internal energy and derivatives wrt pressure and temperature at the qps
      Real e, de_dp, de_dT;
      _fp.e_from_p_T(MetaPhysicL::raw_value(_porepressure[_qp][_phase_num]) * _pressure_to_Pascals,
                     MetaPhysicL::raw_value(_temperature[_qp]) + _t_c2k,
                     e,
                     de_dp,
                     de_dT);
      (*_internal_energy)[_qp] = e;
      (*_dinternal_energy_dp)[_qp] = de_dp * _pressure_to_Pascals;
      (*_dinternal_energy_dT)[_qp] = de_dT;
    }
  }

  if (_compute_enthalpy)
  {
    if (is_ad)
      (*_enthalpy)[_qp] = _fp.h_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                         _temperature[_qp] + _t_c2k);
    else
    {
      // Enthalpy and derivatives wrt pressure and temperature at the qps
      Real h, dh_dp, dh_dT;
      _fp.h_from_p_T(MetaPhysicL::raw_value(_porepressure[_qp][_phase_num]) * _pressure_to_Pascals,
                     MetaPhysicL::raw_value(_temperature[_qp]) + _t_c2k,
                     h,
                     dh_dp,
                     dh_dT);
      (*_enthalpy)[_qp] = h;
      (*_denthalpy_dp)[_qp] = dh_dp * _pressure_to_Pascals;
      (*_denthalpy_dT)[_qp] = dh_dT;
    }
  }
}

template class PorousFlowConcentrationDensityTempl<false>;
template class PorousFlowConcentrationDensityTempl<true>;
