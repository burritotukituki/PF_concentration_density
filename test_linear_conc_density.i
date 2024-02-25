[Mesh]
    [mesh]
      type = GeneratedMeshGenerator
      dim = 2
      nx = 3
      xmin = 0
      xmax = 1
      ny = 100
      ymin = -1
      ymax = 0
      bias_y = 0.925
    []
  []
  
  [GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 -10 0'
  []
  
  [Variables]
    [pp]
      initial_condition = 1E6
    []
    [xnacl]
       initial_condition = 0.05
    []
  []
  
  # [ICs]
  #   [pressure]
  #     type = FunctionIC
  #     function = "-10*1000*y"
  #     variable = pp
  #   []
  #  []
  
   [BCs]
  #  [constant_injection_of_xnacl]
  #    type = PorousFlowSink
  #    variable = xnacl
  #    flux_function = -1E-5
  #    boundary = top
  #  []
    [constant_injection_of_xnacl]
      type = FunctionDirichletBC
      variable = xnacl
      function = 'if(t>=0, 0.3, 0.05)'
      boundary = top
    []
  []

   [AuxVariables]
    [density_xnacl]
       order = CONSTANT
       family = MONOMIAL
    []
   []
    [AuxKernels]
    [density_mass]
       type = PorousFlowPropertyAux
       variable = density_xnacl
       property = density
    []
   []

  [Kernels]
    [mass0]
      type = PorousFlowMassTimeDerivative
      variable = pp
      fluid_component = 0
    []
    [mass1]
      type = PorousFlowMassTimeDerivative
      variable = xnacl
      fluid_component = 1
    []
    [flux0]
      type = PorousFlowAdvectiveFlux
      fluid_component = 0
      variable = pp
    []
    [flux1]
      type = PorousFlowAdvectiveFlux
      fluid_component = 1
      variable = xnacl
    []
  []
  
  [UserObjects]
    [dictator]
      type = PorousFlowDictator
      porous_flow_vars = 'pp xnacl'
      number_fluid_phases = 1
      number_fluid_components = 2
    []
    [pc]
      type = PorousFlowCapillaryPressureConst
      pc = 0
    []
  []
  
  [FluidProperties]
    [brine]
      type = BrineFluidProperties
    []
  []
  
  [Materials]
    [temperature]
      type = PorousFlowTemperature
      temperature = 293
    []
    [ppss]
      type = PorousFlow1PhaseP
      porepressure = pp
      capillary_pressure = pc
    []
    [massfrac]
      type = PorousFlowMassFraction
      mass_fraction_vars = 'xnacl'
    []
    [brine]
      type = PorousFlowBrine
      phase = 0
      xnacl = xnacl
    []
    [permeability]
      type = PorousFlowPermeabilityConst
      permeability = '1e-11 0 0 0 1e-11 0 0 0 3e-14'
    []
    [porosity]
      type = PorousFlowPorosityConst
      porosity = 0.1
    []
    [relperm]
      type = PorousFlowRelativePermeabilityConst
      phase = 0
    []
  []
  
  [Preconditioning]
    active = basic
    [mumps_is_best_for_parallel_jobs]
      type = SMP
      full = true
      petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
      petsc_options_value = ' lu       mumps'
    []
    [basic]
      type = SMP
      full = true
      petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
      petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
    []
    [moose]
      type = SMP
      full = true
      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
      petsc_options_value = 'hypre boomeramg 500'
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = NEWTON
    end_time = 1e5
    start_time = -1E3
    nl_max_its = 25
    l_max_its = 100
    dtmax = 1e4
    nl_abs_tol = 1e-9
    automatic_scaling = true
    compute_scaling_once = false
    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 1E3
      growth_factor = 2
      cutback_factor = 0.5
    []
  []
  
  [Outputs]
    exodus = true
  []



