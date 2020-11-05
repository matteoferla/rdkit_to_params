########################################################################################################################

__doc__ = \
"""
``_PoserMixin`` is a base that adds pyrosetta functionality if avaliable. Adds:

*the bound method ``.test(outfile=None)``, which returns a pose.
*the static method ``.params_to_pose(paramsfile:str, name:str)``, which returns a pose.

"""
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "5 November 2020 A.D."
__license__ = "MIT"
__version__ = "1.1.3"
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
        self.log.debug(f'Testing params {"(writing out as "+outfile+")" if outfile is not None else "(no saving)"}')
        pose = self.to_pose(relax=True)
        assert pose.total_residue() == 1, 'Residue failed.'
        if outfile:
            pose.dump_pdb(outfile)
        return pose

    def to_pose(self, relax=False):
        pose = pyrosetta.Pose()
        # add paramsblock
        rts = self.add_residuetype(pose)
        # add new residue
        lig = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(rts.name_map(self.NAME))
        pose.append_residue_by_jump(lig, 1)
        if relax:
            cycles = 15
            scorefxn = pyrosetta.get_fa_scorefxn()
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
            relax.apply(pose)
        return pose

    def add_residuetype(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.core.chemical.ResidueTypeSet:
        rts = pose.conformation().modifiable_residue_type_set_for_conf(pyrosetta.rosetta.core.chemical.FULL_ATOM_t)
        buffer = pyrosetta.rosetta.std.stringbuf(self.dumps())
        stream = pyrosetta.rosetta.std.istream(buffer)
        new = pyrosetta.rosetta.core.chemical.read_topology_file(stream, self.NAME,
                                                                 rts)  # no idea what the second argument does
        rts.add_base_residue_type(new)
        return rts


    @staticmethod
    def params_to_pose(paramsfile:str, name:str) -> pyrosetta.Pose:
        """
        Staticmethod to get a pose from a params file.

        :param paramsfile: params file.
        :param name: 3-letter residue name.
        :return:
        """
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([paramsfile])
        resiset = pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        lig = pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue( resiset.name_map( name ) )
        pose.append_residue_by_jump(lig, 1)
        return pose
