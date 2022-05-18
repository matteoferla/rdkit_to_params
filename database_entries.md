## Database entries

This module is able to read params and generate a Params instance 
or convert a params file directly to a pose.

```python
pose = Params.load(params_path, skip_unknown=True).to_pose()  # parsed
pose = Params.params_to_pose(params_path, 'LIG')  # no parsing
```
Do note that of the 1520 params files in `pyrosetta/database/chemical/residue_type_sets/fa_standard/residue_types/`
only a fraction can generate a pose successfully:

* 302 without parsing
* 279 with parsing

The difference is mainly due to broken paths.
A couple however segfault.
So the analysis needs Sqlitedict module comes in handy —also used in Fragmenstein.

```python
import os
import pyrosetta
import pyrosetta_help as ph
from rdkit import Chem, RDConfig
from rdkit_to_params import Params
from sqlitedict import SqliteDict
from rdkit.Chem import PandasTools
import pandas as pd

# capture to log
logger = ph.configure_logger()

extra_options= ph.make_option_string(no_optH=False,
                                    ex1=None,
                                    ex2=None,
                                    #mute='all',
                                    ignore_unrecognized_res=False, # raise error!
                                    load_PDB_components=False,
                                    ignore_waters=False)

pyrosetta.init(extra_options=extra_options)

def get_filenames(folder):
    filenames = []
    for filename in os.listdir(folder):
        path = os.path.join(folder, filename)
        if os.path.isdir(path):
            filenames.extend(get_filenames(path))
        elif '.params' in filename:
            filenames.append(path)
    return filenames

folder = os.path.join(os.path.split(pyrosetta.__file__)[0], # noqa
                      'database',
                      'chemical',
                      'residue_type_sets',
                      'fa_standard',
                      'residue_types')  

params_paths = get_filenames(folder)
summary = SqliteDict('aa.db', autocommit=True)

for params_path in params_paths:
    filename = os.path.split(params_path)[1]
    if filename in summary:
        pass
    data = {'error': 'segfault'}
    summary[filename] = data
    try:
        params = Params.load(params_path, skip_unknown=True)
        data['path'] = params_path.replace(folder, '')
        data['name'] = params.NAME
        data['type'] = str(params.TYPE)
        data['comments'] = str(params.comments).replace('#','')
        summary[filename] = data
        pdbblock = ph.get_pdbstr(params.to_pose())
        data['SMILES'] = Chem.MolToSmiles(Chem.MolFromPDBBlock(pdbblock))
        data['error'] = None
    except (MemoryError, RuntimeError) as error:
        data = summary[filename]
        data['error'] = error.args[0].strip().split('\n')[-1]
    except Exception as error:
        data = summary[filename]
        summary[filename]['error'] = f'{error.__class__.__name__}: {error}'
    finally:
        summary[filename] = data
        
from rdkit_to_params import Params

for params_path in params_paths:
    filename = os.path.split(params_path)[1]
    if 'error_original' in summary[filename]:
        continue
    data = summary[filename]
    data['error_original'] = 'segfault'
    data['good_original'] = False
    summary[filename] = data
    try:
        name = data['name']
        pose = Params.params_to_pose(params_path, name)
        data['good_original'] = bool(len(pose.sequence()))
        data['error_original'] = None
    except Exception as error:
        data['error_original'] = error.args[0].strip().split('\n')[-1]
    finally:
        summary[filename] = data

summary_table = pd.DataFrame.from_dict(dict(summary), orient='index')
summary_table['is_name3_too_long'] = summary_table.error.str.contains('is not 3 char long').fillna(False)
summary_table['formatting_error'] = summary_table.error.fillna('')\
                                                 .str.extract(r'(\w+) entry .* is not formatted correctly')
summary_table.to_csv('round_trip.csv')
```

The latter dataset can be found in [examples/round_trip.csv](examples/round_trip.csv)

    ['001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '012', '013', '014', '015', '016', '020', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '202', '203', '204', '205', '207', '208', '209', '210', '211', '2MA', '2MG', '2MU', '302', '303', '304', '305', '306', '307', '313', '314', '315', '316', '317', '318', '320', '332', '333', '3TD', '401', '402', '404', '405', '406', '407', '408', '409', '410', '411', '412', '4OC', '501', '502', '503', '504', '505', '507', '508', '509', '5MC', '601', '621', '623', '631', '633', '6MZ', '701', '703', '704', '7MG', 'A04', 'A06', 'A07', 'A12', 'A20', 'A24', 'A31', 'A32', 'A33', 'A34', 'A43', 'A45', 'A48', 'A68', 'A69', 'A78', 'A80', 'A82', 'A83', 'A84', 'A91', 'A92', 'A94', 'A98', 'ABA', 'AIB', 'ALA', 'ALT', 'APA', 'APN', 'ASN', 'B02', 'B06', 'B12', 'B19', 'B19', 'B21', 'B27', 'B28', 'B30', 'B31', 'B35', 'B36', 'B38', 'B3A', 'B3C', 'B3D', 'B3E', 'B3F', 'B3G', 'B3H', 'B3I', 'B3K', 'B3L', 'B3M', 'B3N', 'B3O', 'B3P', 'B3Q', 'B3R', 'B3S', 'B3T', 'B3V', 'B3W', 'B3X', 'B3Y', 'B40', 'B44', 'B47', 'B48', 'B49', 'B50', 'B54', 'B56', 'B57', 'B58', 'B59', 'B60', 'B61', 'B62', 'B63', 'B67', 'B74', 'B95', 'B96', 'B97', 'BB8', 'BCS', 'BMT', 'BPY', 'BZP', 'C00', 'C01', 'C03', 'C04', 'C05', 'C12', 'C15', 'C16', 'C20', 'C26', 'C27', 'C30', 'C36', 'C41', 'C42', 'C53', 'C54', 'C55', 'C60', 'C61', 'C80', 'C81', 'C83', 'C84', 'C85', 'C86', 'C88', 'C89', 'C89', 'C91', 'C92', 'C93', 'C94', 'C95', 'CAL', 'CO3', 'CPN', 'CYS', 'CYS', 'CYS', 'CYX', 'DAB', 'DPP', 'EMB', 'FE2', 'G7M', 'GDP', 'GLN', 'GLY', 'GNP', 'GPN', 'GTP', 'HHA', 'HIP', 'HIS', 'HLU', 'HOH', 'HP2', 'HPR', 'HYD', 'IGL', 'ILE', 'LEU', 'LMA', 'LYS', 'M2G', 'MEM', 'MET', 'MTP', 'NLU', 'NVL', 'OC3', 'OHH', 'OPH', 'ORN', 'OSi', 'PHA', 'PHE', 'PHO', 'PRO', 'PSU', 'R1A', 'S56', 'SAL', 'SER', 'TES', 'THR', 'TMA', 'TPN', 'TRP', 'TYR', 'UPN', 'UR3', 'V01', 'V02', 'V03', 'V04', 'VAL', 'VDP', 'VHA', 'VLS', 'VOR', 'VSR', 'dhI']

Do note that if PyRosetta is a mock, `os.path.dirname(pyrosetta.__file__)` will raise an exception.