from srwlib import *
from uti_plot import *
import numpy

if not srwl_uti_proc_is_master(): exit()

####################################################
# LIGHT SOURCE

part_beam = SRWLPartBeam()
part_beam.Iavg               = 0.5
part_beam.partStatMom1.x     = -3.3332869187089774e-05
part_beam.partStatMom1.y     = -9.021806076181468e-07
part_beam.partStatMom1.z     = -1.8667500000000001
part_beam.partStatMom1.xp    = 2.5100184866462523e-06
part_beam.partStatMom1.yp    = -7.882119026507647e-07
part_beam.partStatMom1.gamma = 5870.363475323999
part_beam.arStatMom2[0]      = 1.1467411396e-08
part_beam.arStatMom2[1]      = 0.0
part_beam.arStatMom2[2]      = 2.63789068816e-11
part_beam.arStatMom2[3]      = 2.7199979929599998e-11
part_beam.arStatMom2[4]      = 0.0
part_beam.arStatMom2[5]      = 2.3529412449e-12
part_beam.arStatMom2[10]     = 7.920999999999999e-07

magnetic_fields = []
magnetic_fields.append(SRWLMagFldH(1, 'v', 0.187782, 0, 1, 1))
magnetic_structure = SRWLMagFldU(magnetic_fields,0.057,61.5)
magnetic_field_container = SRWLMagFldC([magnetic_structure], array('d', [0]), array('d', [0]), array('d', [0]))

mesh = SRWLRadMesh(_eStart=999.832015419757,
                   _eFin  =999.832015419757,
                   _ne    =1,
                   _xStart=-0.00338,
                   _xFin  =0.00338,
                   _nx    =201,
                   _yStart=-0.01,
                   _yFin  =0.01,
                   _ny    =201,
                   _zStart=34.366)

stk = SRWLStokes()
stk.allocate(1,201,201)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh
wfr.partBeam = part_beam

initial_mesh = deepcopy(wfr.mesh)

####################################################
# BEAMLINE

srw_oe_array = []
srw_pp_array = []
drift_before_oe_0 = SRWLOptD(20.434)
pp_drift_before_oe_0 = [0,0,1.0,1,0,1.2,3.5,1.2,3.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_0)
srw_pp_array.append(pp_drift_before_oe_0)

oe_0=SRWLOptA(_shape='r',
            _ap_or_ob='a',
            _Dx=0.02,
            _Dy=0.017415,
            _x=0.0,
            _y=0.0)
pp_oe_0 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_0)
srw_pp_array.append(pp_oe_0)

drift_before_oe_1 = SRWLOptD(0.2)
pp_drift_before_oe_1 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_1)
srw_pp_array.append(pp_drift_before_oe_1)

oe_1=SRWLOptA(_shape='r',
            _ap_or_ob='a',
            _Dx=0.015,
            _Dy=0.002591,
            _x=0.0,
            _y=0.0)
pp_oe_1 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_1)
srw_pp_array.append(pp_oe_1)

substrate_mirror = SRWLOptMirPl(_size_tang=0.2,
                                _size_sag=0.015,
                                _ap_shape='r',
                                _sim_meth=2,
                                _treat_in_out=1,
                                _nvx=0,
                                _nvy=0.9999160776597537,
                                _nvz=-0.012955216614688831,
                                _tvx=0,
                                _tvy=0.012955216614688831,
                                _x=0.0,
                                _y=0.0)
substrate_mirror.set_dim_sim_meth(_size_tang=0.2,
                                  _size_sag=0.015,
                                  _ap_shape='r',
                                  _sim_meth=2,
                                  _treat_in_out=1)
substrate_mirror.set_orient(_nvx=0,
                            _nvy=0.9999160776597537,
                            _nvz=-0.012955216614688831,
                            _tvx=0,
                            _tvy=0.012955216614688831,
                            _x=0.0,
                            _y=0.0)

oe_2=SRWLOptG(_mirSub=substrate_mirror,
            _m=1,
            _grDen=1800.0,
            _grDen1=0.08997,
            _grDen2=3.004e-06,
            _grDen3=9.73e-11,
            _grDen4=0.0,
            _grAng= 0.0)
pp_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0,0.0,0.08093907,0.99671905,-1.0,0.0]

srw_oe_array.append(oe_2)
srw_pp_array.append(pp_oe_2)

drift_after_oe_2 = SRWLOptD(34.63)
pp_drift_after_oe_2 = [0,0,1.0,2,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_after_oe_2)
srw_pp_array.append(pp_drift_after_oe_2)

oe_3=SRWLOptA(_shape='r',
            _ap_or_ob='a',
            _Dx=0.00916,
            _Dy=0.02,
            _x=0.0,
            _y=0.0)
pp_oe_3 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_3)
srw_pp_array.append(pp_oe_3)

oe_4 = SRWLOptMirEl(_size_tang=0.42,
                  _size_sag=0.02,
                  _p=89.63,
                  _q=8.006,
                  _ang_graz=0.02181661564992912,
                  _ap_shape='r',
                  _sim_meth=2,
                  _treat_in_out=1,
                  _nvx=-0.9997620270799091,
                  _nvy=0,
                  _nvz=-0.02181488503456112,
                  _tvx=0.02181488503456112,
                  _tvy=0,
                  _x=0.0,
                  _y=0.0)
oe_4.set_dim_sim_meth(_size_tang=0.42,
                      _size_sag=0.02,
                      _ap_shape='r',
                      _sim_meth=2,
                      _treat_in_out=1)
oe_4.set_orient(_nvx=-0.9997620270799091,
               _nvy=0,
               _nvz=-0.02181488503456112,
               _tvx=0.02181488503456112,
               _tvy=0,
               _x=0.0,
               _y=0.0)

pp_oe_4 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_4)
srw_pp_array.append(pp_oe_4)

drift_after_oe_4 = SRWLOptD(8.006)
pp_drift_after_oe_4 = [0,0,1.0,4,0,1.0,1.0,3.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_after_oe_4)
srw_pp_array.append(pp_drift_after_oe_4)

oe_5=SRWLOptA(_shape='r',
            _ap_or_ob='a',
            _Dx=0.01,
            _Dy=1e-05,
            _x=0.0,
            _y=0.0)
pp_oe_5 = [0,0,1.0,0,0,0.4,1.0,0.1,1.0,0,0.0,0.0]

srw_oe_array.append(oe_5)
srw_pp_array.append(pp_oe_5)

drift_before_oe_6 = SRWLOptD(6.01)
pp_drift_before_oe_6 = [0,0,1.0,3,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_6)
srw_pp_array.append(pp_drift_before_oe_6)

oe_6=SRWLOptA(_shape='r',
            _ap_or_ob='a',
            _Dx=0.00654,
            _Dy=0.05,
            _x=0.0,
            _y=0.0)
pp_oe_6 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_6)
srw_pp_array.append(pp_oe_6)

oe_7=SRWLOptL(_Fx=0.85072, _Fy=0.85072, _x=0.0, _y=0.0)
pp_oe_7 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_7)
srw_pp_array.append(pp_oe_7)

drift_after_oe_7 = SRWLOptD(0.991)
pp_drift_after_oe_7 = [0,0,1.0,4,0,2.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_after_oe_7)
srw_pp_array.append(pp_drift_after_oe_7)



####################################################
# PROPAGATION

optBL = SRWLOptC(srw_oe_array, srw_pp_array)



####################################################
# MULTI ELECTRON PROPAGATION

radStokesProp = srwl_wfr_emit_prop_multi_e(part_beam,
                                           magnetic_field_container,
                                           initial_mesh,
                                           1,
                                           0.01,
                                           500000,
                                           5,
                                           20,
                                           'output_srw_script_me.dat',
                                           1.0,
                                           optBL,
                                           _char=0)
