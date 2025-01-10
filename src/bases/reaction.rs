use super::{
    atom::Atom,
    initialization::Either,
    molecule::{Molecule, MoleculeBlueprint},
};
use serde::Deserialize;

#[derive(Debug, Clone)]
pub enum ParticleBlueprint {
    Molecule(MoleculeBlueprint),
    Atom(String),
}

#[derive(Debug, Clone)]
pub struct FormingReactionBlueprint {
    pub reactants: Either<ParticleBlueprint>,
    pub products: ParticleBlueprint,
    pub atom_idxs: (Option<usize>, Option<usize>), // atom_idx1 in mol1 and atom_idx2 in mol2 which when they
                                                   // interact they form molecule
}
#[derive(Debug, Clone)]
pub struct BreakingReactionBlueprint {
    pub reactants: ParticleBlueprint,
    pub products: Either<ParticleBlueprint>,
    pub bond_idx: usize,
}
