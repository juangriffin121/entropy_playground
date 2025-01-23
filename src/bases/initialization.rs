use super::{
    molecule::MoleculeBlueprint,
    reaction::{
        BreakingReactionBlueprint, FormingReactionBlueprint, ParticleBlueprint,
        ReactionParticipant, ReactionParticipantInfo,
    },
    recipient::BlueprintAtomIndex,
};
use std::collections::{HashMap, HashSet};

pub fn recursive_definitions(
    particle: &ParticleBlueprint,
    reaction_registry: &mut HashMap<
        MoleculeBlueprint,
        Vec<(FormingReactionBlueprint, BreakingReactionBlueprint)>,
    >,
) {
    match particle {
        ParticleBlueprint::Atom(_name) => {
            // not adding anything to registry
        }
        ParticleBlueprint::Molecule(ref molecule) => {
            if molecule.bonds.is_empty() || molecule.atoms.len() == 1 {
                unreachable!("molecule with no bonds or one atom");
            }
            let new_molecule = if !molecule.canonical {
                molecule.canonicalize().0
            } else {
                molecule.clone()
            };
            let molecule = new_molecule;
            if reaction_registry.contains_key(&molecule) {
                return;
            }
            let mut molecule_reactions = Vec::new();
            for (bond_idx, bond) in molecule.bonds.iter().enumerate() {
                let remaining: Vec<_> = molecule
                    .bonds
                    .clone()
                    .into_iter()
                    .filter(|b| b != bond)
                    .collect();
                let (fragments, breaking_reaction, forming_reaction) =
                    get_reactions_and_fragments(&molecule, remaining, *bond);
                molecule_reactions.push((forming_reaction, breaking_reaction));

                match fragments {
                    SingleOrPair::One(new_particle) => {
                        recursive_definitions(&new_particle, reaction_registry);
                    }
                    SingleOrPair::Two((particle1, particle2)) => {
                        if let ParticleBlueprint::Molecule(ref _molecule) = particle1 {
                            recursive_definitions(&particle1, reaction_registry);
                        }
                        if let ParticleBlueprint::Molecule(ref _molecule) = particle2 {
                            recursive_definitions(&particle2, reaction_registry);
                        }
                    }
                    _ => unreachable!("found more than two fragments {:?}", fragments),
                }
            }
            reaction_registry.insert(molecule.clone(), molecule_reactions);
        }
    }
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub enum Either<T, S> {
    Left(T),
    Right((S)),
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub enum SingleOrPair<T> {
    One(T),
    Two((T, T)),
}

fn get_reactions_and_fragments(
    molecule: &MoleculeBlueprint,
    remaining: Vec<(BlueprintAtomIndex, BlueprintAtomIndex)>,
    bond: (BlueprintAtomIndex, BlueprintAtomIndex),
) -> (
    SingleOrPair<ParticleBlueprint>,
    BreakingReactionBlueprint,
    FormingReactionBlueprint,
) {
    let fragments_info = get_fragments(&molecule.atoms, remaining, bond.0, bond.1);
    let fragments = fragments_info.particles();
    match &fragments_info {
        ReactionParticipantInfo::Single((fragment, new_idx_map)) => {
            let breaking_reaction = BreakingReactionBlueprint {
                reactants: molecule.clone(),
                products_info: fragments_info.clone(),
                broken_bond: bond,
            };

            let forming_reaction = FormingReactionBlueprint {
                products: molecule.clone(),
                reactants_info: ReactionParticipantInfo::Single((
                    fragment.clone(),
                    reverse(&new_idx_map),
                )),
                atom_idxs: (Some(new_idx_map[&bond.0]), Some(new_idx_map[&bond.1])),
            };
            (fragments, breaking_reaction, forming_reaction)
        }
        ReactionParticipantInfo::Pair((participant, participant2)) => {
            let breaking_reaction = BreakingReactionBlueprint {
                reactants: molecule.clone(),
                products_info: fragments_info.clone(),
                broken_bond: bond,
            };

            let (reactant, atom1_id) = match participant {
                ReactionParticipant::Molecule((molecule_blueprint, new_idx_map)) => (
                    ReactionParticipant::Molecule((
                        molecule_blueprint.clone(),
                        reverse(&new_idx_map),
                    )),
                    Some(new_idx_map[&bond.0]),
                ),
                ReactionParticipant::FreeAtom((name, new_atom_pos)) => (
                    ReactionParticipant::FreeAtom((name.clone(), *new_atom_pos)),
                    None,
                ),
            };
            println!("");
            println!("{:#?} {:#?}", reactant, atom1_id);

            let (reactant2, atom2_id) = match participant2 {
                ReactionParticipant::Molecule((molecule_blueprint, new_idx_map)) => (
                    ReactionParticipant::Molecule((
                        molecule_blueprint.clone(),
                        reverse(&new_idx_map),
                    )),
                    Some(new_idx_map[&bond.1]),
                ),
                ReactionParticipant::FreeAtom((name, new_atom_pos)) => (
                    ReactionParticipant::FreeAtom((name.clone(), *new_atom_pos)),
                    None,
                ),
            };
            println!("{:#?} {:#?}", reactant2, atom2_id);
            println!("");

            let forming_reactants_info = ReactionParticipantInfo::Pair((reactant, reactant2));

            let forming_reaction = FormingReactionBlueprint {
                products: molecule.clone(),
                reactants_info: forming_reactants_info,
                atom_idxs: (atom1_id, atom2_id),
            };

            (fragments, breaking_reaction, forming_reaction)
        }
    }
}
fn get_fragments(
    atoms: &Vec<String>,
    remaining: Vec<(BlueprintAtomIndex, BlueprintAtomIndex)>,
    atom1_id: BlueprintAtomIndex,
    atom2_id: BlueprintAtomIndex,
) -> ReactionParticipantInfo {
    let mut atoms_seen: HashSet<BlueprintAtomIndex> = HashSet::new();
    let mut bonds_seen: HashSet<(BlueprintAtomIndex, BlueprintAtomIndex)> = HashSet::new();
    let num_atoms = atoms.len();
    let bonds_map: HashMap<BlueprintAtomIndex, Vec<BlueprintAtomIndex>> = (0..num_atoms)
        .map(|i| {
            (
                BlueprintAtomIndex(i),
                remaining
                    .iter()
                    .filter(|(id1, id2)| i == id1.0 || i == id2.0)
                    .map(|(id1, id2)| if i == id1.0 { *id2 } else { *id1 })
                    .collect::<Vec<BlueprintAtomIndex>>(),
            )
        })
        .collect();

    let current = atom1_id;
    let (atom_idxs, bonds) = dfs(current, &mut atoms_seen, &mut bonds_seen, &bonds_map);
    let participant = make_new_particle(atoms, atom_idxs, bonds);
    if !atoms_seen.contains(&atom2_id) {
        let current = atom2_id;
        let (atom_idxs, bonds) = dfs(current, &mut atoms_seen, &mut bonds_seen, &bonds_map);
        if let Some(_id) = (0..num_atoms).find(|&a| !atoms_seen.contains(&BlueprintAtomIndex(a))) {
            panic!()
        } else {
            let participant2 = make_new_particle(atoms, atom_idxs, bonds);
            return ReactionParticipantInfo::Pair((participant, participant2));
        }
    } else {
        if let Some(_id) = (0..num_atoms).find(|&a| !atoms_seen.contains(&BlueprintAtomIndex(a))) {
            panic!()
        } else {
            if let ReactionParticipant::Molecule((molecule, new_idx_map)) = participant {
                return ReactionParticipantInfo::Single((molecule, new_idx_map));
            } else {
                panic!("atom as only fragment")
            }
        }
    }
}

fn dfs(
    current: BlueprintAtomIndex,
    atoms_seen: &mut HashSet<BlueprintAtomIndex>,
    bonds_seen: &mut HashSet<(BlueprintAtomIndex, BlueprintAtomIndex)>,
    bonds_map: &HashMap<BlueprintAtomIndex, Vec<BlueprintAtomIndex>>,
) -> (
    Vec<BlueprintAtomIndex>,
    Vec<(BlueprintAtomIndex, BlueprintAtomIndex)>,
) {
    let bonds = &bonds_map[&current];
    let mut molecule_bonds = Vec::new();
    let mut molecule_atoms = Vec::new();

    atoms_seen.insert(current);
    molecule_atoms.push(current);

    for &bonded in bonds {
        let bond = (
            BlueprintAtomIndex(current.0.min(bonded.0)),
            BlueprintAtomIndex(current.0.max(bonded.0)),
        );
        if !bonds_seen.contains(&bond) {
            bonds_seen.insert(bond);
            if !atoms_seen.contains(&bonded) {
                let (next_atoms, next_bonds) = dfs(bonded, atoms_seen, bonds_seen, bonds_map);
                molecule_atoms.extend(next_atoms);
                molecule_bonds.extend(next_bonds);
            }
            molecule_bonds.push(bond);
        }
    }

    (molecule_atoms, molecule_bonds)
}

fn make_new_molecule(
    atoms: &Vec<String>,
    atom_idxs: Vec<BlueprintAtomIndex>,
    bonds: Vec<(BlueprintAtomIndex, BlueprintAtomIndex)>,
) -> (
    MoleculeBlueprint,
    HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
) {
    let new_mol_atoms: Vec<String> = atom_idxs.iter().map(|&i| atoms[i.0].clone()).collect();
    let new_idx_map: HashMap<BlueprintAtomIndex, BlueprintAtomIndex> = atom_idxs
        .iter()
        .enumerate()
        .map(|(new_index, &old_index)| (old_index, BlueprintAtomIndex(new_index)))
        .collect();
    let new_bonds: Vec<(BlueprintAtomIndex, BlueprintAtomIndex)> = bonds
        .iter()
        .map(|(i1, i2)| (new_idx_map[i1], new_idx_map[i2]))
        .collect();
    let molecule = MoleculeBlueprint {
        atoms: new_mol_atoms,
        bonds: new_bonds,
        canonical: false,
    };

    let (canonical_molecule, canonical_map) = molecule.canonicalize();

    // Update new_idx_map to reflect changes in indices after canonicalization
    let new_idx_map = new_idx_map
        .iter()
        .map(|(&old_index, &pre_canonical_index)| (old_index, canonical_map[&pre_canonical_index]))
        .collect();
    (canonical_molecule, new_idx_map)
}

fn make_new_particle(
    atoms: &Vec<String>,
    atom_idxs: Vec<BlueprintAtomIndex>,
    bonds: Vec<(BlueprintAtomIndex, BlueprintAtomIndex)>,
) -> ReactionParticipant {
    if atom_idxs.len() > 1 {
        let (molecule, new_idx_map) = make_new_molecule(&atoms, atom_idxs, bonds);
        ReactionParticipant::Molecule((molecule, new_idx_map))
    } else {
        ReactionParticipant::FreeAtom((atoms[atom_idxs[0].0].clone(), atom_idxs[0]))
    }
}

fn reverse(
    idx_map: &HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
) -> HashMap<BlueprintAtomIndex, BlueprintAtomIndex> {
    (idx_map
        .iter()
        .map(|(&i, &j)| (j, i))
        .collect::<HashMap<BlueprintAtomIndex, BlueprintAtomIndex>>())
}
