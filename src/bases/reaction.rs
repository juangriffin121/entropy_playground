use std::collections::HashMap;

use crate::bases::bond::order_bond;

use super::{
    atom::{Atom, Element},
    bond::{add_bond, Bond, BondProperties},
    initialization::Either,
    molecule::{Molecule, MoleculeBlueprint},
    recipient::{self, Recipient},
};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum ParticleBlueprint {
    Molecule(MoleculeBlueprint),
    Atom(String),
}

#[derive(Debug, Clone)]
pub struct FormingReactionBlueprint {
    pub reactants: Either<ParticleBlueprint>,
    pub products: MoleculeBlueprint,
    pub atom_idxs: (Option<usize>, Option<usize>), // atom_idx1 in mol1 and atom_idx2 in mol2 which when they
    // interact they form molecule
    pub new_idx_map: Either<Option<HashMap<usize, usize>>>,
}

#[derive(Debug, Clone)]
pub struct BreakingReactionBlueprint {
    pub reactants: MoleculeBlueprint,
    pub products: Either<ParticleBlueprint>,
    pub broken_bond: (usize, usize),
    pub new_idx_map: Either<Option<HashMap<usize, usize>>>,
}

pub fn check_forming_reactions<'a>(
    atom1: &Atom,
    atom2: &Atom,
    forming_reactions: &'a Option<
        HashMap<
            (Either<ParticleBlueprint>, (Option<usize>, Option<usize>)),
            FormingReactionBlueprint,
        >,
    >,
    molecule_registry: &HashMap<usize, Molecule>,
    molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
) -> Option<(Either<usize>, &'a FormingReactionBlueprint)> {
    let possible_reactions = forming_reactions.as_ref()?;
    let (idx1, particle1) = get_particle_from_atom(atom1, molecule_registry, molecule_blueprints);
    let (idx2, particle2) = get_particle_from_atom(atom2, molecule_registry, molecule_blueprints);
    let (idxs, reactants) = match (&particle1, &particle2) {
        (ParticleBlueprint::Molecule(molecule1), ParticleBlueprint::Molecule(molecule2)) => {
            if idx1 == idx2 {
                (Either::One(idx1), Either::One(particle1))
            } else {
                (
                    Either::Two((idx1, idx2)),
                    Either::Two((particle1, particle2)),
                )
            }
        }
        (_, _) => (
            Either::Two((idx1, idx2)),
            Either::Two((particle1, particle2)),
        ),
    };

    let reaction = possible_reactions.get(&(reactants, (atom1.molecule_id, atom2.molecule_id)));
    Some((idxs, reaction?))
}

fn get_particle_from_atom(
    atom: &Atom,
    molecule_registry: &HashMap<usize, Molecule>,
    molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
) -> (usize, ParticleBlueprint) {
    match atom.molecule_id {
        Some(idx) => (
            idx,
            ParticleBlueprint::Molecule({
                let blueprint_id = molecule_registry
                    .get(&idx)
                    .expect("atom with refernce to molecule that doesnt exists")
                    .blueprint_id;

                molecule_blueprints
                    .as_ref()
                    .expect("reactions but no molecules")[blueprint_id]
                    .clone()
            }),
        ),
        None => (
            atom.id,
            ParticleBlueprint::Atom(atom.chemical_nature.clone()),
        ),
    }
}

impl FormingReactionBlueprint {
    pub fn apply(
        &self,
        atom_ids: (usize, usize),
        reactant_idxs: Either<usize>,
        atoms: &mut Vec<Atom>,
        molecule_registry: &mut HashMap<usize, Molecule>,
        bond_map: &mut HashMap<(usize, usize), Bond>,
        next_molecule_id: &mut usize,
        blueprint_ids_map: &HashMap<MoleculeBlueprint, usize>,
        bond_blueprints: &HashMap<(String, String), BondProperties>,
    ) {
        let ordered = order_bond(atom_ids);

        let (left, right) = atoms.split_at_mut(ordered.1);
        let a = &mut left[ordered.0];
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

        let properties = bond_blueprints
            .get(&name_pair)
            .expect("chemical pair isnt in bond_blueprints");
        let (k, equilibrium_distance, breaking_distance) = (
            properties.strength,
            properties.equilibrium_distance,
            properties.breaking_distance,
        );
        add_bond(bond_map, a, b, k, equilibrium_distance, breaking_distance);
        match reactant_idxs {
            Either::One(reactant_id) => {
                let mut new_bonds = vec![ordered];
                if let Some(reactant) = molecule_registry.remove(&reactant_id) {
                    new_bonds.extend(reactant.bonds.clone());
                } else {
                    unreachable!(
                        "forming reactions with only one reactant can only happen if reactant is a molecule, no molecule with id {} in registry",
                        reactant_id
                        )
                }

                let blueprint_id = blueprint_ids_map[&self.products];
                let new_mol = Molecule {
                    id: *next_molecule_id,
                    bonds: new_bonds,
                    blueprint_id,
                };
                molecule_registry.insert(*next_molecule_id, new_mol);
                *next_molecule_id += 1;
            }
            Either::Two((reactant_id1, reactant_id2)) => {
                let mut new_bonds = vec![ordered];
                if let Some(reactant) = molecule_registry.remove(&reactant_id1) {
                    new_bonds.extend(reactant.bonds.clone());
                }

                if let Some(reactant) = molecule_registry.remove(&reactant_id1) {
                    new_bonds.extend(reactant.bonds.clone());
                }

                let blueprint_id = blueprint_ids_map[&self.products];
                let new_mol = Molecule {
                    id: *next_molecule_id,
                    bonds: new_bonds,
                    blueprint_id,
                };
                molecule_registry.insert(*next_molecule_id, new_mol);
                *next_molecule_id += 1;
            }
        }
    }
}

pub fn check_breaking_reactions<'a>(
    atom1: &Atom,
    atom2: &Atom,
    breaking_reactions: &'a Option<
        HashMap<(MoleculeBlueprint, (usize, usize)), BreakingReactionBlueprint>,
    >,
    molecule_registry: &HashMap<usize, Molecule>,
    molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
) -> Option<&'a BreakingReactionBlueprint> {
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
    let blueprint = &blueprints[molecule.blueprint_id];

    // HORRIBLE CLONE AND DEREF ILL CHANGE THIS AT SOME POINT
    let reaction = possible_reactions.get(&(blueprint.clone(), (atom1_id_in_mol, atom2_id_in_mol)));
    reaction
}

impl BreakingReactionBlueprint {
    pub fn apply(
        &self,
        atom1_id: usize,
        atom2_id: usize,
        atoms: &mut Vec<Atom>,
        molecule_registry: &mut HashMap<usize, Molecule>,
        bond_map: &mut HashMap<(usize, usize), Bond>,
        next_molecule_id: &mut usize,
        blueprint_ids_map: &HashMap<MoleculeBlueprint, usize>,
    ) {
    }
}
