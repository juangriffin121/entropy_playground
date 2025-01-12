use std::collections::HashMap;

use super::{
    atom::{Atom, Element},
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
}

#[derive(Debug, Clone)]
pub struct BreakingReactionBlueprint {
    pub reactants: MoleculeBlueprint,
    pub products: Either<ParticleBlueprint>,
    pub bond_idx: usize,
}

pub fn check_reactions<'a>(
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
    ) {
    }
}
