./pythonsh prepare_receptor4.py -r 1GX4_atom.clean.pdb -A 'hydrogens'
./pythonsh prepare_ligand4.py -l 1GX4_hetatom.pdb -Z
./pythonsh prepare_gpf4.py -l 1GX4_hetatom.pdbqt -r 1GX4_atom.clean.pdbqt -i 1GX4_atom.clean.pdb_conf.txt
./pythonsh prepare_dpf4.py -l 1GX4_hetatom.pdbqt -r 1GX4_atom.clean.pdbqt
./autogrid4 -p 1GX4_atom.clean.gpf -l 1GX4_hetatom.pdb.glg
./autodock4 -p 1GX4_hetatom_1GX4_atom.dpf -l 1GX4_atom_dock.dlg
echo "1 Over"
