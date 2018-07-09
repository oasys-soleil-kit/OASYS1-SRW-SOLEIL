import sys

from PyQt5.QtWidgets import QApplication
from orangewidget import gui
from orangewidget.settings import Setting
from oasys.widgets import gui as oasysgui
from oasys.widgets import congruence

from soleil.wofrysrw.storage_ring.light_sources.srw_infrared_light_source import SRWIRBendingMagnetLightSource
from soleil.wofrysrw.storage_ring.magnetic_structures.srw_infrared_bending_magnet import SRWIRBendingMagnet

from orangecontrib.srw.widgets.gui.ow_srw_source import OWSRWSource

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet

class OWSRWIRBendingMagnet(OWSRWSource):

    name = "IR Bending Magnet"
    description = "SRW Source: IR Bending Magnet"
    icon = "icons/bending_magnet.png"
    priority = 1

    magnetic_radius = Setting(5.8164)
    magnetic_field = Setting(1.72)
    length = Setting(0.525)

    center_of_straight_section = Setting(6.2025)
    transition_steepness = Setting(60)
    z_start = Setting(-8.0)
    z_end = Setting(8.0)

    want_main_area=1

    def __init__(self):
        super().__init__()

        left_box_2 = oasysgui.widgetBox(self.tab_source, "ID Parameters", addSpace=True, orientation="vertical", height=200)

        oasysgui.lineEdit(left_box_2, self, "magnetic_radius", "Magnetic Radius [m]", labelWidth=260, valueType=float, orientation="horizontal", callback=self.calculateMagneticField)
        oasysgui.lineEdit(left_box_2, self, "magnetic_field", "Magnetic Field [T]", labelWidth=260, valueType=float, orientation="horizontal", callback=self.calculateMagneticRadius)
        oasysgui.lineEdit(left_box_2, self, "length", "Length [m]", labelWidth=260, valueType=float, orientation="horizontal")

        oasysgui.lineEdit(left_box_2, self, "center_of_straight_section", "Center of Straight Section [m]", labelWidth=260, valueType=float, orientation="horizontal")
        oasysgui.lineEdit(left_box_2, self, "transition_steepness", "Transition Steepness", labelWidth=260, valueType=float, orientation="horizontal")

        box = oasysgui.widgetBox(left_box_2, "", orientation="horizontal")

        oasysgui.lineEdit(box, self, "z_start", "By: z start [m]", valueType=float, orientation="horizontal")
        oasysgui.lineEdit(box, self, "z_end", "end [m]", valueType=float, orientation="horizontal")

        gui.rubber(self.controlArea)
        gui.rubber(self.mainArea)

    def get_default_initial_z(self):
        return 0.0 # initial Longitudinal Coordinate

    def get_srw_source(self, electron_beam):

        self.magnetic_radius = BendingMagnet.calculate_magnetic_radius(self.magnetic_field, electron_beam.electron_energy_in_GeV) if self.magnetic_radius == 0.0 else self.magnetic_radius
        self.magnetic_field = BendingMagnet.calculate_magnetic_field(self.magnetic_radius, electron_beam.electron_energy_in_GeV) if self.magnetic_field == 0.0 else self.magnetic_field

        return SRWIRBendingMagnetLightSource(electron_beam=electron_beam,
                                             bending_magnet_magnetic_structure=SRWIRBendingMagnet(radius=self.magnetic_radius,
                                                                                                  magnetic_field=self.magnetic_field,
                                                                                                  length=self.length,
                                                                                                  center_of_straight_section=self.center_of_straight_section,
                                                                                                  transition_steepness=self.transition_steepness,
                                                                                                  z_start=self.z_start,
                                                                                                  z_end=self.z_end,
                                                                                                  n_points=self.wf_number_of_points_for_trajectory_calculation))

    def print_specific_infos(self, srw_source):
        pass

    def get_automatic_sr_method(self):
        return 2

    def get_source_length(self):
        return self.length

    def checkLightSourceSpecificFields(self):
        congruence.checkStrictlyPositiveNumber(self.magnetic_radius, "Magnetic Radius")
        congruence.checkStrictlyPositiveNumber(self.magnetic_field, "Magnetic Field")
        congruence.checkStrictlyPositiveNumber(self.length, "Length")
        congruence.checkPositiveNumber(self.transition_steepness, "Transition Steepness")
        congruence.checkGreaterThan(self.z_end, self.z_start, "z end", "z start")

    def calculateMagneticField(self):
        if self.magnetic_radius > 0:
           self.magnetic_field=BendingMagnet.calculate_magnetic_field(self.magnetic_radius,
                                                                      self.electron_energy_in_GeV)

    def calculateMagneticRadius(self):
        if self.magnetic_field > 0:
           self.magnetic_radius=BendingMagnet.calculate_magnetic_radius(self.magnetic_field,
                                                                        self.electron_energy_in_GeV)

    def receive_specific_syned_data(self, data):
        if isinstance(data._light_source._magnetic_structure, BendingMagnet):
            light_source = data._light_source

            self.magnetic_field = light_source._magnetic_structure._magnetic_field
            self.magnetic_radius = light_source._magnetic_structure._radius
            self.length = light_source._magnetic_structure._length
        else:
            raise ValueError("Syned data not correct")


if __name__ == "__main__":
    a = QApplication(sys.argv)
    ow = OWSRWIRBendingMagnet()
    ow.show()
    a.exec_()
    ow.saveSettings()
