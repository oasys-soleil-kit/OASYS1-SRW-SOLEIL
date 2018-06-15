from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam

from soleil.wofrysrw.storage_ring.magnetic_structures.srw_infrared_bending_magnet import SRWIRBendingMagnet

class SRWIRBendingMagnetLightSource(SRWLightSource):

    def __init__(self,
                 name="Undefined",
                 electron_beam=SRWElectronBeam(),
                 bending_magnet_magnetic_structure=SRWIRBendingMagnet()):


        super().__init__(name,
                         electron_beam=electron_beam,
                         magnetic_structure=bending_magnet_magnetic_structure)
