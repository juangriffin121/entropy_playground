use super::{
    molecule::MoleculeBlueprint,
    reaction::{BreakingReactionBlueprint, FormingReactionBlueprint, ParticleBlueprint},
};
use std::{
    collections::{HashMap, HashSet},
    usize,
};

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
                molecule.canonicalize()
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
                    Either::One(new_particle) => {
                        recursive_definitions(&new_particle, reaction_registry);
                    }
                    Either::Two((particle1, particle2)) => {
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
pub enum Either<T> {
    One(T),
    Two((T, T)),
}

fn get_reactions_and_fragments(
    molecule: &MoleculeBlueprint,
    remaining: Vec<(usize, usize)>,
    bond: (usize, usize),
) -> (
    Either<ParticleBlueprint>,
    BreakingReactionBlueprint,
    FormingReactionBlueprint,
) {
    match get_fragments(&molecule.atoms, remaining, bond.0, bond.1) {
        Either::One((fragment, Some(new_idx_map))) => {
            let breaking_reaction = BreakingReactionBlueprint {
                reactants: molecule.clone(),
                products: Either::One(fragment.clone()),
                broken_bond: bond,
                new_idx_map: Either::One(Some(new_idx_map.clone())),
            };

            let forming_reaction = FormingReactionBlueprint {
                products: molecule.clone(),
                reactants: Either::One(fragment.clone()),
                atom_idxs: (Some(new_idx_map[&bond.0]), Some(new_idx_map[&bond.1])),
                new_idx_map: Either::One(reverse(&Some(new_idx_map))),
            };
            let fragments = Either::One(fragment);
            (fragments, breaking_reaction, forming_reaction)
        }
        Either::One((_, None)) => {
            unreachable!("if breaking a molecule returns only one particle it cant have no bonds")
        }
        Either::Two(((fragment, fragment_idx_map), (fragment2, fragment2_idx_map))) => {
            let breaking_reaction = BreakingReactionBlueprint {
                reactants: molecule.clone(),
                products: Either::Two((fragment.clone(), fragment2.clone())),
                broken_bond: bond,
                new_idx_map: Either::Two((fragment_idx_map.clone(), fragment2_idx_map.clone())),
            };

            let forming_reaction = FormingReactionBlueprint {
                products: molecule.clone(),
                reactants: Either::Two((fragment.clone(), fragment2.clone())),
                atom_idxs: (
                    match fragment_idx_map {
                        Some(ref map) => Some(
                            *map.get(&bond.0)
                                .expect("if reactant1 is a molecule it should contain atom1"),
                        ),
                        None => None,
                    },
                    match fragment2_idx_map {
                        Some(ref map) => Some(
                            *map.get(&bond.1)
                                .expect("if reactant2 is a molecule it should contain atom2"),
                        ),
                        None => None,
                    },
                ),
                new_idx_map: Either::Two((reverse(&fragment_idx_map), reverse(&fragment2_idx_map))),
            };
            let fragments = Either::Two((fragment, fragment2));
            (fragments, breaking_reaction, forming_reaction)
        }
    }
}
fn get_fragments(
    atoms: &Vec<String>,
    remaining: Vec<(usize, usize)>,
    atom1_id: usize,
    atom2_id: usize,
) -> Either<(ParticleBlueprint, Option<HashMap<usize, usize>>)> {
    let mut atoms_seen: HashSet<usize> = HashSet::new();
    let mut bonds_seen: HashSet<(usize, usize)> = HashSet::new();
    let num_atoms = atoms.len();
    let bonds_map: HashMap<usize, Vec<usize>> = (0..num_atoms)
        .map(|i| {
            (
                i,
                remaining
                    .iter()
                    .filter(|(id1, id2)| i == *id1 || i == *id2)
                    .map(|(id1, id2)| if i == *id1 { *id2 } else { *id1 })
                    .collect::<Vec<usize>>(),
            )
        })
        .collect();

    let current = atom1_id;
    let (atom_idxs, bonds) = dfs(current, &mut atoms_seen, &mut bonds_seen, &bonds_map);
    let (particle, particle_idx_map) = make_new_particle(atoms, atom_idxs, bonds);
    if !atoms_seen.contains(&atom2_id) {
        let current = atom2_id;
        let (atom_idxs, bonds) = dfs(current, &mut atoms_seen, &mut bonds_seen, &bonds_map);
        if let Some(_id) = (0..num_atoms).find(|a| !atoms_seen.contains(a)) {
            panic!()
        } else {
            return Either::Two((
                (particle, particle_idx_map),
                make_new_particle(atoms, atom_idxs, bonds),
            ));
        }
    } else {
        if let Some(_id) = (0..num_atoms).find(|a| !atoms_seen.contains(a)) {
            panic!()
        } else {
            return Either::One((particle, particle_idx_map));
        }
    }
}

fn dfs(
    current: usize,
    atoms_seen: &mut HashSet<usize>,
    bonds_seen: &mut HashSet<(usize, usize)>,
    bonds_map: &HashMap<usize, Vec<usize>>,
) -> (Vec<usize>, Vec<(usize, usize)>) {
    let bonds = &bonds_map[&current];
    let mut molecule_bonds = Vec::new();
    let mut molecule_atoms = Vec::new();

    atoms_seen.insert(current);
    molecule_atoms.push(current);

    for &bonded in bonds {
        let bond = (current.min(bonded), current.max(bonded));
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
    atom_idxs: Vec<usize>,
    bonds: Vec<(usize, usize)>,
) -> (MoleculeBlueprint, HashMap<usize, usize>) {
    let new_mol_atoms: Vec<String> = atom_idxs.iter().map(|&i| atoms[i].clone()).collect();
    let new_idx_map: HashMap<usize, usize> = atom_idxs
        .iter()
        .enumerate()
        .map(|(new_index, &old_index)| (old_index, new_index))
        .collect();
    let new_bonds: Vec<(usize, usize)> = bonds
        .iter()
        .map(|(i1, i2)| (new_idx_map[i1], new_idx_map[i2]))
        .collect();
    let molecule = MoleculeBlueprint {
        atoms: new_mol_atoms,
        bonds: new_bonds,
        canonical: false,
    };
    (molecule.canonicalize(), new_idx_map)
}

fn make_new_particle(
    atoms: &Vec<String>,
    atom_idxs: Vec<usize>,
    bonds: Vec<(usize, usize)>,
) -> (ParticleBlueprint, Option<HashMap<usize, usize>>) {
    if atom_idxs.len() > 1 {
        let (molecule, new_idx_map) = make_new_molecule(&atoms, atom_idxs, bonds);
        (ParticleBlueprint::Molecule(molecule), Some(new_idx_map))
    } else {
        (ParticleBlueprint::Atom(atoms[atom_idxs[0]].clone()), None)
    }
}

fn reverse(idx_map: &Option<HashMap<usize, usize>>) -> Option<HashMap<usize, usize>> {
    Some(idx_map.as_ref()?.iter().map(|(&i, &j)| (j, i)).collect())
}
