use super::{
    atom::{Atom, Element},
    bond::{add_bond, remove_bond, Bond, BondProperties},
    initialization::{Either, SingleOrPair},
    molecule::{self, Molecule, MoleculeBlueprint},
    recipient::{BlueprintAtomIndex, MoleculeAtomIndex, MoleculeId, Recipient, WorldAtomId},
};
use crate::bases::bond::order_bond;
use std::{collections::HashMap, process::id};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum ParticleBlueprint {
    Molecule(MoleculeBlueprint),
    Atom(String),
}

#[derive(Debug, Clone)]
pub struct FormingReactionBlueprint {
    pub reactants_info: ReactionParticipantInfo,
    pub products: MoleculeBlueprint,
    pub atom_idxs: (Option<BlueprintAtomIndex>, Option<BlueprintAtomIndex>), // atom_idx1 in mol1 and atom_idx2 in mol2 which when they
                                                                             // interact they form molecule
}

#[derive(Debug, Clone)]
pub enum ReactionParticipant {
    Molecule(
        (
            MoleculeBlueprint,
            HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
        ),
    ),
    FreeAtom((String, BlueprintAtomIndex)),
}

#[derive(Debug, Clone)]
pub enum ReactionParticipantInfo {
    Single(
        (
            MoleculeBlueprint,
            HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
        ),
    ),
    Pair(((ReactionParticipant, ReactionParticipant))),
}

impl FormingReactionBlueprint {
    pub fn key(&self) -> FormingReactionKey {
        let idxs = self.atom_idxs;
        let info = self.reactants_info.clone();
        match (info, idxs) {
            (ReactionParticipantInfo::Pair((p1, p2)), (id1, id2)) => {
                let p1 = match (p1, id1) {
                    (ReactionParticipant::Molecule((mol, _)), Some(id1)) => {
                        ReactionParticipantKey::Molecule((mol, id1))
                    }
                    (ReactionParticipant::FreeAtom((name, _)), None) => {
                        ReactionParticipantKey::FreeAtom(name)
                    }
                    _ => unreachable!("molecule without id or atom with"),
                };
                let p2 = match (p2, id2) {
                    (ReactionParticipant::Molecule((mol, _)), Some(id2)) => {
                        ReactionParticipantKey::Molecule((mol, id2))
                    }
                    (ReactionParticipant::FreeAtom((name, _)), None) => {
                        ReactionParticipantKey::FreeAtom(name)
                    }
                    _ => unreachable!("molecule without id or atom with"),
                };
                FormingReactionKey::Pair((p1, p2))
            }

            (ReactionParticipantInfo::Single((mol, _)), (Some(id1), Some(id2))) => {
                FormingReactionKey::Single((mol, (id1, id2)))
            }
            _ => unreachable!("single molecule without some of the idxs"),
        }
    }
}

impl ReactionParticipantInfo {
    pub fn particles(&self) -> SingleOrPair<ParticleBlueprint> {
        match self {
            Self::Single((molecule, _)) => {
                SingleOrPair::One(ParticleBlueprint::Molecule(molecule.clone()))
            }
            Self::Pair((participant1, participant2)) => {
                let p1 = match participant1 {
                    ReactionParticipant::Molecule((molecule, _)) => {
                        ParticleBlueprint::Molecule(molecule.clone())
                    }
                    ReactionParticipant::FreeAtom((atom, _)) => {
                        ParticleBlueprint::Atom(atom.clone())
                    }
                };

                let p2 = match participant2 {
                    ReactionParticipant::Molecule((molecule, _)) => {
                        ParticleBlueprint::Molecule(molecule.clone())
                    }
                    ReactionParticipant::FreeAtom((atom, _)) => {
                        ParticleBlueprint::Atom(atom.clone())
                    }
                };
                SingleOrPair::Two((p1, p2))
            }
        }
    }
}
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum ReactionParticipantKey {
    Molecule((MoleculeBlueprint, BlueprintAtomIndex)),
    FreeAtom(String),
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum FormingReactionKey {
    Single((MoleculeBlueprint, (BlueprintAtomIndex, BlueprintAtomIndex))),
    Pair(((ReactionParticipantKey, ReactionParticipantKey))),
}

#[derive(Debug, Clone)]
pub struct BreakingReactionBlueprint {
    pub reactants: MoleculeBlueprint,
    pub products_info: ReactionParticipantInfo,
    pub broken_bond: (BlueprintAtomIndex, BlueprintAtomIndex),
}

pub fn check_forming_reactions<'a>(
    atom1: &Atom,
    atom2: &Atom,
    forming_reactions: &'a Option<HashMap<FormingReactionKey, FormingReactionBlueprint>>,
    molecule_registry: &HashMap<MoleculeId, Molecule>,
    molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
) -> Option<(
    Either<
        MoleculeId,
        (
            Either<MoleculeId, WorldAtomId>,
            Either<MoleculeId, WorldAtomId>,
        ),
    >,
    &'a FormingReactionBlueprint,
)> {
    let possible_reactions = forming_reactions.as_ref()?;
    let (idx1, particle1) = get_particle_from_atom(atom1, molecule_registry, molecule_blueprints);
    let (idx2, particle2) = get_particle_from_atom(atom2, molecule_registry, molecule_blueprints);
    let key = match (
        (&particle1, atom1.id_in_molecule),
        (&particle2, atom2.id_in_molecule),
    ) {
        (
            (ParticleBlueprint::Molecule(molecule1), Some(id_in_mol1)),
            (ParticleBlueprint::Molecule(molecule2), Some(id_in_mol2)),
        ) => {
            if idx1 == idx2 {
                FormingReactionKey::Single((
                    molecule1.clone(),
                    (id_in_mol1.to_bp(), id_in_mol2.to_bp()),
                ))
            } else {
                FormingReactionKey::Pair((
                    ReactionParticipantKey::Molecule((molecule1.clone(), id_in_mol1.to_bp())),
                    ReactionParticipantKey::Molecule((molecule2.clone(), id_in_mol2.to_bp())),
                ))
            }
        }
        (
            (ParticleBlueprint::Molecule(molecule1), Some(id_in_mol1)),
            (ParticleBlueprint::Atom(name), None),
        ) => FormingReactionKey::Pair((
            ReactionParticipantKey::Molecule((molecule1.clone(), id_in_mol1.to_bp())),
            ReactionParticipantKey::FreeAtom(name.clone()),
        )),
        (
            (ParticleBlueprint::Atom(name), None),
            (ParticleBlueprint::Molecule(molecule2), Some(id_in_mol2)),
        ) => FormingReactionKey::Pair((
            ReactionParticipantKey::FreeAtom(name.clone()),
            ReactionParticipantKey::Molecule((molecule2.clone(), id_in_mol2.to_bp())),
        )),
        ((ParticleBlueprint::Atom(name1), None), (ParticleBlueprint::Atom(name2), None)) => {
            FormingReactionKey::Pair((
                ReactionParticipantKey::FreeAtom(name1.clone()),
                ReactionParticipantKey::FreeAtom(name2.clone()),
            ))
        }
        _ => unreachable!("atom with mol_id or molecule without"),
    };

    let reaction = possible_reactions.get(&key);
    // println!("key: {:#?}\nreaction: {:#?}\n", key, reaction);
    let idxs: Either<
        MoleculeId,
        (
            Either<MoleculeId, WorldAtomId>,
            Either<MoleculeId, WorldAtomId>,
        ),
    > = match (idx1, idx2) {
        (Either::Left(mol_id1), Either::Left(mol_id2)) => {
            if mol_id1 == mol_id2 {
                Either::Left(mol_id1)
            } else {
                Either::Right((Either::Left(mol_id1), Either::Left(mol_id2)))
            }
        }
        (Either::Left(mol_id1), Either::Right(atom_id2)) => {
            Either::Right((Either::Left(mol_id1), Either::Right(atom_id2)))
        }
        (Either::Right(atom_id1), Either::Left(mol_id2)) => {
            Either::Right((Either::Right(atom_id1), Either::Left(mol_id2)))
        }
        (Either::Right(atom_id1), Either::Right(atom_id2)) => {
            Either::Right((Either::Right(atom_id1), Either::Right(atom_id2)))
        }
    };
    Some((idxs, reaction?))
}

fn get_particle_from_atom(
    atom: &Atom,
    molecule_registry: &HashMap<MoleculeId, Molecule>,
    molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
) -> (Either<MoleculeId, WorldAtomId>, ParticleBlueprint) {
    match atom.molecule_id {
        Some(idx) => (
            Either::Left(idx),
            ParticleBlueprint::Molecule({
                let blueprint_id = molecule_registry
                    .get(&idx)
                    .expect(&format!(
                        "atom with refernce to molecule that doesnt exists {:#?}\n{:#?}",
                        atom, molecule_registry
                    ))
                    .blueprint_id;

                molecule_blueprints
                    .as_ref()
                    .expect("reactions but no molecules")
                    .get(blueprint_id)
                    .expect("no blueprint with that id")
                    .clone()
            }),
        ),
        None => (
            Either::Right(atom.id),
            ParticleBlueprint::Atom(atom.chemical_nature.clone()),
        ),
    }
}

// might not be changing atom1_id_in_mol
impl FormingReactionBlueprint {
    pub fn apply(
        &self,
        atom_ids: (WorldAtomId, WorldAtomId),
        reactant_idxs: Either<
            MoleculeId,
            (
                Either<MoleculeId, WorldAtomId>,
                Either<MoleculeId, WorldAtomId>,
            ),
        >,
        atoms: &mut Vec<Atom>,
        molecule_registry: &mut HashMap<MoleculeId, Molecule>,
        bond_map: &mut HashMap<(WorldAtomId, WorldAtomId), Bond>,
        next_molecule_id: &mut MoleculeId,
        blueprint_ids_map: &HashMap<MoleculeBlueprint, usize>,
        bond_blueprints: &HashMap<(String, String), BondProperties>,
    ) {
        let ordered_atom_ids = order_bond(atom_ids);

        let (left, right) = atoms.split_at_mut(ordered_atom_ids.1 .0);
        let a = &mut left[ordered_atom_ids.0 .0];
        let b = &mut right[0];
        let name_pair = {
            let a_name = a.chemical_nature.clone();
            let b_name = b.chemical_nature.clone();
            if a_name > b_name {
                (b_name, a_name)
            } else {
                (a_name, b_name)
            }
        };

        let properties = bond_blueprints.get(&name_pair).expect(&format!(
            "chemical pair isnt in bond_blueprints {:?} ",
            name_pair,
        ));
        let (k, equilibrium_distance, breaking_distance) = (
            properties.strength,
            properties.equilibrium_distance,
            properties.breaking_distance,
        );
        add_bond(bond_map, a, b, k, equilibrium_distance, breaking_distance);
        let atom1 = &atoms[atom_ids.0 .0];
        let atom2 = &atoms[atom_ids.1 .0];
        match (reactant_idxs, &self.reactants_info) {
            (
                Either::Left(molecule_id),
                ReactionParticipantInfo::Single((molecule_blueprint, new_idx_map)),
            ) => {
                let atom1_pos = atom1
                    .id_in_molecule
                    .expect("atom1 not in molecule when in branch of single reactor");
                let atom2_pos = atom2
                    .id_in_molecule
                    .expect("atom2 not in molecule when in branch of single reactor");

                let molecule = molecule_registry
                    .remove(&molecule_id)
                    .expect("no molecule with that id in registry");
                let new_idx_map =
                    check_and_map(new_idx_map, &molecule, &molecule_blueprint, &atoms);

                let molecule_atoms = molecule.atoms.clone();

                let atom1_new_pos = new_idx_map
                    .get(&atom1_pos)
                    .expect("atom1_pos not in new_idx_map");
                let atom2_new_pos = new_idx_map
                    .get(&atom1_pos)
                    .expect("atom1_pos not in new_idx_map");
                let mut new_bonds = vec![(*atom1_new_pos, *atom2_new_pos)];

                let (world_atom_ids, new_atoms_pos, molecule_bonds) =
                    process_molecule_reactant(&new_idx_map, molecule);

                new_bonds.extend(molecule_bonds);
                let world_atom_ids = order_array(&world_atom_ids, &new_atoms_pos);

                let blueprint_id = blueprint_ids_map
                    .get(&self.products)
                    .expect("product_blueprint not in blueprint_ids_map");
                let new_mol = Molecule {
                    id: *next_molecule_id,
                    atoms: world_atom_ids,
                    bonds_: new_bonds,
                    blueprint_id: *blueprint_id,
                };
                molecule_registry.insert(*next_molecule_id, new_mol);
                update_atom_molecule_refs2(
                    atoms,
                    *next_molecule_id,
                    &new_atoms_pos,
                    &molecule_atoms,
                );

                next_molecule_id.0 += 1;
            }
            (
                Either::Right((reactant_id1, reactant_id2)),
                ReactionParticipantInfo::Pair((reactant1_info, reactant2_info)),
            ) => {
                let (mut par1_atoms, mut new_atoms_pos1, par1_bonds, atom1_new_pos) =
                    process_particle_reactant(
                        reactant1_info,
                        reactant_id1,
                        molecule_registry,
                        atoms,
                        atom1,
                    );

                let (par2_atoms, new_atoms_pos2, par2_bonds, atom2_new_pos) =
                    process_particle_reactant(
                        reactant2_info,
                        reactant_id2,
                        molecule_registry,
                        atoms,
                        atom2,
                    );
                update_atom_molecule_refs2(atoms, *next_molecule_id, &new_atoms_pos1, &par1_atoms);
                update_atom_molecule_refs2(atoms, *next_molecule_id, &new_atoms_pos2, &par2_atoms);

                let mut new_bonds = vec![(atom1_new_pos, atom2_new_pos)];
                new_bonds.extend(par1_bonds);
                new_bonds.extend(par2_bonds);
                par1_atoms.extend(par2_atoms);
                new_atoms_pos1.extend(new_atoms_pos2);

                let world_atom_ids = order_array(&par1_atoms, &new_atoms_pos1);

                let blueprint_id = blueprint_ids_map
                    .get(&self.products)
                    .expect("product_blueprint not in blueprint_ids_map 2nd case");
                let new_mol = Molecule {
                    id: *next_molecule_id,
                    atoms: world_atom_ids,
                    bonds_: new_bonds,
                    blueprint_id: *blueprint_id,
                };
                if !check_molecule(&new_mol, &self.products, atoms) {
                    panic!("molecule didnt pass check {:?}", new_mol)
                }
                molecule_registry.insert(*next_molecule_id, new_mol);

                next_molecule_id.0 += 1;
            }
            (Either::Left(_), ReactionParticipantInfo::Pair(_)) => {
                unreachable!("single reactant idx but pair of reactant info")
            }
            (Either::Right(_), ReactionParticipantInfo::Single(_)) => {
                unreachable!("single reactant but found two particles")
            }
        }
    }
}

pub fn check_breaking_reactions<'a>(
    atom1: &Atom,
    atom2: &Atom,
    breaking_reactions: &'a Option<
        HashMap<
            (MoleculeBlueprint, (BlueprintAtomIndex, BlueprintAtomIndex)),
            BreakingReactionBlueprint,
        >,
    >,
    molecule_registry: &HashMap<MoleculeId, Molecule>,
    molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
) -> Option<(MoleculeId, &'a BreakingReactionBlueprint)> {
    /// Need to check molecule is valid with its blueprint
    let possible_reactions = breaking_reactions.as_ref()?;
    let blueprints = molecule_blueprints.as_ref()?;
    if atom1.molecule_id != atom2.molecule_id {
        unreachable!("atoms in bond that just broke arent in the same molecule")
    }
    let molecule_id = atom1
        .molecule_id
        .expect("atom with bond that broke not in molecule");
    let molecule = molecule_registry
        .get(&molecule_id)
        .expect("molecule_id in atoms that broke a bond isnt in the registry");

    let atom1_id_in_mol = atom1
        .id_in_molecule
        .expect("atom in bond that broke doesnt have id_in_molecule");
    let atom2_id_in_mol = atom2
        .id_in_molecule
        .expect("atom in bond that broke doesnt have id_in_molecule");
    let blueprint = blueprints
        .get(molecule.blueprint_id)
        .expect("product_blueprint not in blueprint_ids_map 2nd case");

    /// HORRIBLE CLONE AND DEREF ILL CHANGE THIS AT SOME POINT
    let reaction = possible_reactions.get(&(
        blueprint.clone(),
        (
            BlueprintAtomIndex(atom1_id_in_mol.0),
            BlueprintAtomIndex(atom2_id_in_mol.0),
        ),
    ));
    Some((molecule_id, reaction?))
}

impl BreakingReactionBlueprint {
    pub fn apply(
        &self,
        atom1_id: WorldAtomId,
        atom2_id: WorldAtomId,
        atoms: &mut Vec<Atom>,
        molecule_id: MoleculeId,
        molecule_registry: &mut HashMap<MoleculeId, Molecule>,
        bond_map: &mut HashMap<(WorldAtomId, WorldAtomId), Bond>,
        next_molecule_id: &mut MoleculeId,
        molecule_blueprints: &Vec<MoleculeBlueprint>,
        blueprint_ids_map: &HashMap<MoleculeBlueprint, usize>,
    ) {
        let ordered_atom_ids = order_bond((atom1_id, atom2_id));
        let (left, right) = atoms.split_at_mut(ordered_atom_ids.1 .0);
        let a = &mut left[ordered_atom_ids.0 .0];
        let b = &mut right[0];
        remove_bond(bond_map, a, b);

        let molecule = molecule_registry
            .remove(&molecule_id)
            .expect("molecule_id in atoms that broke a bond isnt in the registry");
        let blueprint = &molecule_blueprints[molecule.blueprint_id];
        if !check_molecule(&molecule, blueprint, atoms) {
            panic!("molecule found not valid")
        };

        match &self.products_info {
            ReactionParticipantInfo::Single((product_blueprint, position_mapping)) => {
                let position_mapping = transform_position_mapping(position_mapping);
                let new_atom_positions = molecule
                    .atoms
                    .iter()
                    .enumerate()
                    .map(|(old_atom_position, _atom_id)| {
                        *position_mapping
                            .get(&MoleculeAtomIndex(old_atom_position))
                            .expect("position_mapping doesnt have old_atom_position")
                    })
                    .collect();
                let new_atoms_vec = order_array(&molecule.atoms, &new_atom_positions);
                let blueprint_id = blueprint_ids_map[product_blueprint];
                let new_bonds: Vec<(MoleculeAtomIndex, MoleculeAtomIndex)> = molecule
                    .bonds_
                    .iter()
                    .map(|(old_id1, old_id2)| {
                        (
                            *position_mapping
                                .get(old_id1)
                                .expect("position_mapping doesnt have old_id1"),
                            *position_mapping
                                .get(old_id2)
                                .expect("position_mapping doesnt have old_id2"),
                        )
                    })
                    .collect();

                molecule_registry.insert(
                    *next_molecule_id,
                    Molecule {
                        id: *next_molecule_id,
                        bonds_: new_bonds,
                        blueprint_id,
                        atoms: new_atoms_vec,
                    },
                );

                update_atom_molecule_refs2(
                    atoms,
                    *next_molecule_id,
                    &new_atom_positions,
                    &molecule.atoms,
                );

                next_molecule_id.0 += 1;
            }
            ReactionParticipantInfo::Pair((product1_info, product2_info)) => {
                match product1_info {
                    ReactionParticipant::Molecule((product_blueprint, position_mapping)) => {
                        let position_mapping = transform_position_mapping(position_mapping);
                        let new_atom_positions = molecule
                            .atoms
                            .iter()
                            .enumerate()
                            .filter(|(old_atom_position, _)| {
                                position_mapping
                                    .contains_key(&MoleculeAtomIndex(*old_atom_position))
                            })
                            .map(|(old_atom_position, _atom_id)| {
                                position_mapping[&MoleculeAtomIndex(old_atom_position)]
                            })
                            .collect();
                        let molecule_atoms: Vec<WorldAtomId> = molecule
                            .atoms
                            .iter()
                            .enumerate()
                            .filter(|(old_atom_position, _)| {
                                position_mapping
                                    .contains_key(&MoleculeAtomIndex(*old_atom_position))
                            })
                            .map(|(_old_atom_position, atom_id)| atom_id.clone())
                            .collect();
                        let new_atoms_vec = order_array(&molecule.atoms, &new_atom_positions);
                        let blueprint_id = blueprint_ids_map
                            .get(product_blueprint)
                            .expect("blueprint_ids_map doesnt have product_blueprint");
                        let new_bonds: Vec<(MoleculeAtomIndex, MoleculeAtomIndex)> = molecule
                            .bonds_
                            .iter()
                            .filter(|(old_id1, old_id2)| {
                                position_mapping.contains_key(old_id1)
                                    && position_mapping.contains_key(old_id2)
                            })
                            .map(|(old_id1, old_id2)| {
                                (position_mapping[old_id1], position_mapping[old_id2])
                            })
                            .collect();

                        molecule_registry.insert(
                            *next_molecule_id,
                            Molecule {
                                id: *next_molecule_id,
                                bonds_: new_bonds,
                                blueprint_id: *blueprint_id,
                                atoms: new_atoms_vec,
                            },
                        );
                        update_atom_molecule_refs2(
                            atoms,
                            *next_molecule_id,
                            &new_atom_positions,
                            &molecule_atoms,
                        );

                        next_molecule_id.0 += 1;
                    }
                    ReactionParticipant::FreeAtom((nature, new_pos)) => {
                        //check bonds is empty
                        //check nature matches atom nature
                        let mut_atom1 = &mut atoms[atom1_id.0];
                        mut_atom1.molecule_id = None;
                        mut_atom1.id_in_molecule = None;
                    }
                };

                match product2_info {
                    ReactionParticipant::Molecule((product_blueprint, position_mapping)) => {
                        let position_mapping = transform_position_mapping(position_mapping);
                        let new_atom_positions = molecule
                            .atoms
                            .iter()
                            .enumerate()
                            .filter(|(old_atom_position, _)| {
                                position_mapping
                                    .contains_key(&MoleculeAtomIndex(*old_atom_position))
                            })
                            .map(|(old_atom_position, _atom_id)| {
                                position_mapping[&MoleculeAtomIndex(old_atom_position)]
                            })
                            .collect();
                        let molecule_atoms: Vec<WorldAtomId> = molecule
                            .atoms
                            .iter()
                            .enumerate()
                            .filter(|(old_atom_position, _)| {
                                position_mapping
                                    .contains_key(&MoleculeAtomIndex(*old_atom_position))
                            })
                            .map(|(_old_atom_position, atom_id)| atom_id.clone())
                            .collect();
                        let new_atoms_vec = order_array(&molecule_atoms, &new_atom_positions);
                        let blueprint_id = blueprint_ids_map
                            .get(product_blueprint)
                            .expect("blueprint_ids_map doesnt have product_blueprint");
                        let new_bonds: Vec<(MoleculeAtomIndex, MoleculeAtomIndex)> = molecule
                            .bonds_
                            .iter()
                            .filter(|(old_id1, old_id2)| {
                                position_mapping.contains_key(old_id1)
                                    && position_mapping.contains_key(old_id2)
                            })
                            .map(|(old_id1, old_id2)| {
                                (position_mapping[old_id1], position_mapping[old_id2])
                            })
                            .collect();

                        molecule_registry.insert(
                            *next_molecule_id,
                            Molecule {
                                id: *next_molecule_id,
                                bonds_: new_bonds,
                                blueprint_id: *blueprint_id,
                                atoms: new_atoms_vec,
                            },
                        );
                        update_atom_molecule_refs2(
                            atoms,
                            *next_molecule_id,
                            &new_atom_positions,
                            &molecule_atoms,
                        );

                        next_molecule_id.0 += 1;
                    }
                    ReactionParticipant::FreeAtom((nature, new_pos)) => {
                        let mut_atom2 = &mut atoms[atom2_id.0];
                        mut_atom2.molecule_id = None;
                        mut_atom2.id_in_molecule = None;
                    }
                }
            }
        }
    }
}

fn process_particle_reactant(
    reactant_blueprint: &ReactionParticipant,
    reactant_id: Either<MoleculeId, WorldAtomId>,
    molecule_registry: &mut HashMap<MoleculeId, Molecule>,
    atoms: &Vec<Atom>,
    atom_in_bond: &Atom,
) -> (
    Vec<WorldAtomId>,
    Vec<MoleculeAtomIndex>,
    Vec<(MoleculeAtomIndex, MoleculeAtomIndex)>,
    MoleculeAtomIndex,
) {
    // println!(
    //     "reactant_bp: {:?}\n reactant_id: {:?}\n atom_in_bond_id: {:?}",
    //     reactant_blueprint, reactant_id, atom_in_bond.id_in_molecule
    // );
    match (reactant_blueprint, reactant_id, atom_in_bond.id_in_molecule) {
        (
            ReactionParticipant::Molecule((molecule_blueprint, new_idx_map)),
            Either::Left(molecule_id),
            Some(atom_in_bond_pos),
        ) => {
            let molecule = molecule_registry
                .remove(&molecule_id)
                .expect("no molecule with that id in registry");
            let new_idx_map = check_and_map(new_idx_map, &molecule, &molecule_blueprint, &atoms);
            let (atoms, new_atoms_pos, new_bonds) =
                process_molecule_reactant(&new_idx_map, molecule);
            let atom_in_bond_new_pos = new_idx_map[&atom_in_bond_pos];

            (atoms, new_atoms_pos, new_bonds, atom_in_bond_new_pos)
        }
        (ReactionParticipant::FreeAtom((_name, new_atom_pos)), Either::Right(atom_id), None) => (
            vec![atom_id],
            vec![MoleculeAtomIndex(new_atom_pos.0)],
            Vec::new(),
            MoleculeAtomIndex(new_atom_pos.0),
        ),
        _ => unreachable!("incompatible types of new_idx_map, particle or reactant_id"),
    }
}

fn process_molecule_reactant(
    new_idx_map: &HashMap<MoleculeAtomIndex, MoleculeAtomIndex>,
    molecule: Molecule,
) -> (
    Vec<WorldAtomId>,
    Vec<MoleculeAtomIndex>,
    Vec<(MoleculeAtomIndex, MoleculeAtomIndex)>,
) {
    let molecule_bonds = molecule
        .bonds_
        .iter()
        .map(|(old_id1, old_id2)| (new_idx_map[old_id1], new_idx_map[old_id2]))
        .collect();
    let new_atoms_pos = molecule
        .atoms
        .iter()
        .enumerate()
        .map(|(old_atom_id, _atom)| new_idx_map[&MoleculeAtomIndex(old_atom_id)])
        .collect();

    (molecule.atoms, new_atoms_pos, molecule_bonds)
}
fn order_array(
    array: &Vec<WorldAtomId>,
    new_positions: &Vec<MoleculeAtomIndex>,
) -> Vec<WorldAtomId> {
    if array.len() != new_positions.len() {
        panic!(
            "different array sizes to order_array {:#?} {:#?}",
            array, new_positions
        );
    }

    let mut shuffled = vec![array[0].clone(); array.len()]; // Create a new array of the same length
    for (i, &new_pos) in new_positions.iter().enumerate() {
        shuffled[new_pos.0] = array[i].clone();
    }
    shuffled
}

fn check_molecule(molecule: &Molecule, blueprint: &MoleculeBlueprint, atoms: &Vec<Atom>) -> bool {
    /// CHECK BONDS ARE THE SAME TOO
    molecule
        .atoms
        .iter()
        .zip(&blueprint.atoms)
        .all(|(atom_id, atom_nature)| atoms[atom_id.0].chemical_nature == *atom_nature)
}

fn transform_position_mapping(
    new_idx_map: &HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
) -> HashMap<MoleculeAtomIndex, MoleculeAtomIndex> {
    new_idx_map
        .iter()
        .map(|(old_id, new_id)| (MoleculeAtomIndex(old_id.0), MoleculeAtomIndex(new_id.0)))
        .collect()
}

fn check_and_map(
    new_idx_map: &HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
    molecule: &Molecule,
    blueprint: &MoleculeBlueprint,
    atoms: &Vec<Atom>,
) -> HashMap<MoleculeAtomIndex, MoleculeAtomIndex> {
    let molecule_is_valid = check_molecule(molecule, blueprint, atoms);
    if molecule_is_valid {
        transform_position_mapping(new_idx_map)
    } else {
        panic!("invalid molecule")
    }
}

fn update_atom_molecule_refs(
    atoms: &mut [Atom],
    bonding_atom_ids: &(WorldAtomId, WorldAtomId),
    ordered_ids: (WorldAtomId, WorldAtomId),
    next_molecule_id: MoleculeId,
    new_positions: (MoleculeAtomIndex, MoleculeAtomIndex),
) {
    let (left, right) = atoms.split_at_mut(ordered_ids.1 .0);
    let a = &mut left[ordered_ids.0 .0];
    let b = &mut right[0];
    a.molecule_id = Some(next_molecule_id);
    b.molecule_id = Some(next_molecule_id);

    let mut_atoms = {
        if a.id == bonding_atom_ids.0 {
            (a, b)
        } else {
            (b, a)
        }
    };
    mut_atoms.0.id_in_molecule = Some(new_positions.0);
    mut_atoms.1.id_in_molecule = Some(new_positions.1);
}

fn update_atom_molecule_refs2(
    atoms: &mut [Atom],
    next_molecule_id: MoleculeId,
    new_atom_positions: &Vec<MoleculeAtomIndex>,
    molecule_atoms: &Vec<WorldAtomId>,
) {
    for (atom_id, new_id) in molecule_atoms.iter().zip(new_atom_positions) {
        atoms[atom_id.0].molecule_id = Some(next_molecule_id);
        atoms[atom_id.0].id_in_molecule = Some(*new_id);
    }
}
