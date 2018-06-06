import numpy
from srwlib import SRWLMagFldM, SRWLMagFld3D, array

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet

from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructure

class SRWIRBendingMagnet(BendingMagnet, SRWMagneticStructure):

    def __init__(self,
				 radius = 0.0,
                 magnetic_field = 1.72,
				 center_of_straight_section = 12.405/2,
				 length = 1.05/2,
				 transition_steepness = 60,
				 z_start = -8.0,
				 z_end = 8.0,
				 n_points = 16001,
                 ):

        BendingMagnet.__init__(self, radius, magnetic_field, length)

        self.center_of_straight_section = center_of_straight_section
        self.transition_steepness = transition_steepness
        self.z_start = z_start
        self.z_end = z_end
        self.n_points = n_points

    def get_SRWMagneticStructure(self):
        z_points = numpy.linspace(self.z_start, self.z_end, self.n_points)

        B_x = numpy.zeros(self.n_points)
        B_z = numpy.zeros(self.n_points)
        B_y = self._magnetic_field * (((((z_points + self.center_of_straight_section)**self.transition_steepness)/(self._length**self.transition_steepness + (z_points + self.center_of_straight_section)**self.transition_steepness))-1) +
									  ((((z_points - self.center_of_straight_section)**self.transition_steepness)/(self._length**self.transition_steepness + (z_points - self.center_of_straight_section)**self.transition_steepness))-1))

        return SRWLMagFld3D(_arBx=array('d', B_x), _arBy=array('d', B_y), _arBz=array('d', B_z),
							_nx=1,                 _ny=1,                 _nz=self.n_points,
							_rx=0.0,               _ry=0.0,               _rz=numpy.abs(self.z_end-self.z_start),
							_nRep=1)

    def to_python_code_aux(self):
        text_code = "magnetic_structure = SRWLMagFldM(" + str(self._magnetic_field) + ", 1, 'n', " + str(self._length) + ")"

        return text_code
