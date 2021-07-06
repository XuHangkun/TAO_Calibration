"""
    Chi2 of anchor optimizer
    ~~~~~~~~~~~~~~~~~~~~~~

    :author: Xu Hangkun (许杭锟)
    :copyright: © 2021 Xu Hangkun <xuhangkun@ihep.ac.cn>
    :license: MIT, see LICENSE for more details.
"""

from utils.create_new_map import create_map_with_ideal_calib_line

class AnchorOptimizeChi2:
    """Chi2 which can be used to optimize anchor position
    """

    def __init__(self,ideal_map,symmetry=True,max_radius=650):
        self.ideal_map = ideal_map
        self.symmetry = symmetry
        self.max_radius = max_radius

    def __call__(self,theta_1,theta_2,phi_2):
        new_map = create_map_with_ideal_calib_line(
            self.ideal_map,{"theta":theta_1,"phi":0},
            {"theta":theta_2,"phi":phi_2},
            symmetry = self.symmetry
        )
        chi2 = self.ideal_map.diff_chi2(new_map,
            max_radius=self.max_radius)
        return chi2
