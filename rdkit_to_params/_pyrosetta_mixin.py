########################################################################################################################

__doc__ = \
"""
``_PoserMixin`` is a base that adds pyrosetta functionality if avaliable. Adds:

*the bound method ``.test(outfile=None)``, which returns a pose.
*the static method ``.params_to_pose(paramsfile:str, name:str)``, which returns a pose.

"""
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "25 June 2020 A.D."
__license__ = "MIT"
__version__ = "1.0.4"
__citation__ = "None."

########################################################################################################################


import pyrosetta, os

from typing import Optional

class _PoserMixin:
    def test(self, outfile:Optional[str]=None) -> pyrosetta.Pose:
        """
        Makes a pose with the ligand to see if it works.

        :param outfile: optionally save file.
        :return: pose
        """

        self.log.debug(f'Testing params! {outfile}')
        self.dump('_test.params')
        pose = self.__class__.params_to_pose('_test.params', self.NAME)
        os.remove('_test.params')
        if outfile:
            pose.dump_pdb(outfile)
        return pose

    @staticmethod
    def params_to_pose(paramsfile:str, name:str) -> pyrosetta.Pose:
        """
        Staticmethod to get a pose from a params file.

        :param paramsfile: params file.
        :param name: 3-letter residue name.
        :return:
        """
        pose = pyrosetta.rosetta.core.pose.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([paramsfile])
        resiset = pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        ## This is the most convoluted way of getting it. But I don't known any otherway.
        v = pyrosetta.rosetta.core.chemical.ResidueTypeFinder(resiset).get_all_possible_residue_types()
        ligtype = [vv for vv in v if vv.name3() == name][0]
        lig = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(ligtype)
        pose.append_residue_by_jump(lig, 1)
        cycles = 15
        scorefxn = pyrosetta.get_fa_scorefxn()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.apply(pose)
        return pose
