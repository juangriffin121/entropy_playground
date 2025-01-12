use super::{
    atom::{Atom, DistributionJSON, Element, ElementsJSON},
    bond::{Bond, BondForce, BondProperties},
    checks::JSONChecker,
    initialization::Either,
    loader::load,
    molecule::{Molecule, MoleculeBlueprint, MoleculesJSON},
    physics::{Force, Gravity, Grid, GridJSON},
};
use crate::bases::{
    bond::BondDefinitions,
    initialization::recursive_definitions,
    reaction::{BreakingReactionBlueprint, FormingReactionBlueprint, ParticleBlueprint},
};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::{
    collections::{HashMap, HashSet},
    usize,
};
use std::{fs, hash::Hash};

pub struct Recipient {
    pub shape: (f32, f32), //can change it to Shape trait to implement shapes other than rectangle
    pub iterations: i32,
    pub dt: f32,
    pub simulation_rate: i32,
    pub grid: Option<Grid>,

    pub elements_data: HashMap<String, Element>,
    pub contents: Vec<Atom>,

    pub molecule_registry: HashMap<usize, Molecule>,
    pub molecule_blueprints: Option<Vec<MoleculeBlueprint>>,
    pub next_molecule_id: usize,

    pub bond_registry: HashMap<usize, Bond>,
    pub bond_blueprints: Option<HashMap<(String, String), BondProperties>>,
    pub next_bond_id: usize,

    pub forming_reactions: Option<
        HashMap<
            (Either<ParticleBlueprint>, (Option<usize>, Option<usize>)),
            FormingReactionBlueprint,
        >,
    >,
    pub breaking_reactions: Option<HashMap<(MoleculeBlueprint, usize), BreakingReactionBlueprint>>,

    pub force_registry: Vec<Box<dyn Force>>,
}

#[derive(Clone, serde::Serialize, serde::Deserialize, Debug)]
pub struct RecipientJSON {
    pub shape: (f32, f32),
    pub position_file: String,
    pub energy: Option<f32>,
    pub iterations: i32,
    pub dt: f32,
    pub simulation_rate: i32,
    pub references: References,
    pub use_grid: bool,
}

#[derive(Clone, serde::Serialize, serde::Deserialize, Debug)]
pub struct References {
    elements: String,
    distribution: Option<String>,
    molecule_blueprints: Option<String>,
    bond_blueprints: Option<String>,
}

impl Recipient {
    pub fn from_json(filepath: &str, rng: &mut impl Rng) -> Recipient {
        let data = RecipientJSON::from_file(filepath).expect("Unable to load RecipientJSON");
        let (elements_data, distribution, molecule_blueprints, bond_blueprints) =
            data.references.load_and_check();

        let (forming_reactions, breaking_reactions, molecule_blueprints) =
            match all_fragmentations(molecule_blueprints) {
                Some((a, b, c)) => (Some(a), Some(b), Some(c)),
                None => (None, None, None),
            };

        let shape = data.shape;
        let paths = load(&data.position_file).expect("couldnt load paths from position_file");

        let num_atoms: usize = paths.iter().map(|(_, points)| points.len()).sum();
        let mut atoms = Vec::new();
        let mut id_: usize = 0;

        let biggest_radius_f32 = elements_data
            .values()
            .max_by(|x, y| {
                x.radius
                    .partial_cmp(&y.radius)
                    .expect("atoms with None radius")
            })
            .map(|x| x.radius)
            .expect("couldnt find biggest_radius");

        let biggest_radius: u32 = biggest_radius_f32.ceil() as u32;
        for (i, path) in paths.iter().enumerate() {
            let (_closed, points) = path;
            let points = remove_duplicates(points.clone(), biggest_radius);

            for point in points {
                let chemical_nature = match distribution {
                    Some(ref dist) => nature_from_distribution(&dist, rng),
                    None => nature_without_distribution(i),
                };
                atoms.push(Atom::new(
                    id_.clone(),
                    (point.0 as f32, point.1 as f32),
                    initialize_velocity(data.energy.unwrap() / num_atoms as f32, 1.0),
                    &chemical_nature,
                    &elements_data,
                ));
                id_ += 1;
            }
        }

        let grid = if data.use_grid {
            let cell_size = (2.0 * biggest_radius_f32, 2.0 * biggest_radius_f32);
            let cells = (
                (shape.0 / cell_size.0).ceil() as i32,
                (shape.1 / cell_size.1).ceil() as i32,
            );
            Some(Grid::new(cells, cell_size))
        } else {
            None
        };
        Recipient {
            shape,
            iterations: data.iterations,
            dt: data.dt,
            simulation_rate: data.simulation_rate,
            grid,

            contents: atoms,
            elements_data,

            molecule_registry: HashMap::new(),
            molecule_blueprints,
            next_molecule_id: 0,

            bond_registry: HashMap::new(),
            bond_blueprints,
            next_bond_id: 0,

            forming_reactions,
            breaking_reactions,
            force_registry: vec![Box::new(Gravity { g: 1.0 }), Box::new(BondForce {})],
        }
    }

    pub fn add_bond(
        &mut self,
        atom1_id: usize,
        atom2_id: usize,
        k: f32,
        equilibrium_distance: f32,
        breaking_distance: f32,
    ) {
        let bond = Bond::new(
            atom1_id,
            atom2_id,
            k,
            equilibrium_distance,
            breaking_distance,
        );
        self.bond_registry.insert(self.next_bond_id, bond);
        self.contents[atom1_id].bonds.insert(self.next_bond_id);
        self.contents[atom2_id].bonds.insert(self.next_bond_id);
        self.next_bond_id += 1;
    }

    pub fn remove_bond(&mut self, bond_id: usize) {
        if let Some(bond) = self.bond_registry.get(&bond_id) {
            let atom1_id = bond.atom1_id;
            self.contents[atom1_id].bonds.remove(&bond_id);
            let atom2_id = bond.atom2_id;
            self.contents[atom2_id].bonds.remove(&bond_id);

            for molecule in &mut self.molecule_registry {
                unimplemented!()
                // deal with later
            }
        }

        self.bond_registry.remove(&bond_id);
    }
}

impl RecipientJSON {
    pub fn from_file(path: &str) -> serde_json::Result<RecipientJSON> {
        let data = fs::read_to_string(path).expect("Unable to read file");
        let rec: RecipientJSON = serde_json::from_str(&data)?;
        Ok(rec)
    }
}
impl References {
    fn load_and_check(
        &self,
    ) -> (
        HashMap<String, Element>,
        Option<Vec<(String, f32)>>,
        Option<HashMap<String, MoleculeBlueprint>>,
        Option<HashMap<(String, String), BondProperties>>,
    ) {
        let elements_data = ElementsJSON::load(&self.elements);
        let distribution = match &self.distribution {
            Some(path) => Some(DistributionJSON::load(path)),
            None => None,
        };
        let molecule_blueprints = match &self.molecule_blueprints {
            Some(path) => Some(MoleculesJSON::load(path)),
            None => None,
        };
        let bond_blueprints = match &self.bond_blueprints {
            Some(path) => Some(BondDefinitions::load(path)),
            None => None,
        };

        JSONChecker::all_checks(
            &elements_data,
            &distribution,
            &molecule_blueprints,
            &bond_blueprints,
        );

        let elements_data = elements_data.data;

        let distribution = match distribution {
            Some(d) => Some(d.to_cdf()),
            None => None,
        };
        let molecule_blueprints = match molecule_blueprints {
            Some(m) => Some(m.blueprints),
            None => None,
        };
        let bond_blueprints = match bond_blueprints {
            Some(b) => Some(b.to_tuple_map()),
            None => None,
        };

        (
            elements_data,
            distribution,
            molecule_blueprints,
            bond_blueprints,
        )
    }
}
fn remove_duplicates(points: Vec<(f64, f64)>, biggest_radius: u32) -> Vec<(f32, f32)> {
    let scale = 1.0 / (biggest_radius as f32);
    let unique_points: HashSet<_> = points
        .into_iter()
        .map(|(x, y)| {
            (
                (x * scale as f64).round() as i32,
                (y * scale as f64).round() as i32,
            )
        })
        .collect();
    unique_points
        .into_iter()
        .map(|(x, y)| (x as f32 / scale as f32, y as f32 / scale as f32))
        .collect()
}

fn initialize_velocity(energy_per_atom: f32, mass: f32) -> (f32, f32) {
    let velocity_magnitude = (2.0 * energy_per_atom / mass).sqrt();

    let normal_dist =
        Normal::new(0.0, velocity_magnitude as f64).expect("Failed creating normal distribution");
    let mut rng = rand::thread_rng();

    let vx = normal_dist.sample(&mut rng) as f32;
    let vy = normal_dist.sample(&mut rng) as f32;

    (vx, vy)
}

fn nature_from_distribution(distribution: &Vec<(String, f32)>, rng: &mut impl Rng) -> String {
    let random_val: f32 = rng.gen();
    for (element, cumulative_prob) in distribution {
        if random_val <= *cumulative_prob {
            return element.to_string();
        }
    }
    unreachable!()
}

fn nature_without_distribution(i: usize) -> String {
    let chemical_nature = match i % 4 {
        0 => "C",
        1 => "N",
        2 => "O",
        3 => "H",
        _ => unreachable!("Modulo 4 should only produce 0, 1, 2, or 3"),
    };
    // let chemical_nature = ChemicalNature::Hidrogen;
    chemical_nature.to_string()
}

fn all_fragmentations(
    explicit_molecules: Option<HashMap<String, MoleculeBlueprint>>,
) -> Option<(
    HashMap<(Either<ParticleBlueprint>, (Option<usize>, Option<usize>)), FormingReactionBlueprint>,
    HashMap<(MoleculeBlueprint, usize), BreakingReactionBlueprint>,
    Vec<MoleculeBlueprint>,
)> {
    let mut reaction_registry = HashMap::new();

    for (molecule_name, molecule) in explicit_molecules.as_ref()?.iter() {
        recursive_definitions(
            &ParticleBlueprint::Molecule(molecule.clone()),
            &mut reaction_registry,
        );
    }

    let forming_reactions: HashMap<
        (Either<ParticleBlueprint>, (Option<usize>, Option<usize>)),
        FormingReactionBlueprint,
    > = reaction_registry
        .iter()
        .map(|(_, molecule_reactions)| {
            molecule_reactions
                .iter()
                .map(|(reaction, _)| {
                    (
                        (reaction.reactants.clone(), reaction.atom_idxs),
                        reaction.clone(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .flatten()
        .collect();
    let breaking_reactions: HashMap<(MoleculeBlueprint, usize), BreakingReactionBlueprint> =
        reaction_registry
            .iter()
            .map(|(_, molecule_reactions)| {
                molecule_reactions
                    .iter()
                    .map(|(_, reaction)| {
                        (
                            (reaction.reactants.clone(), reaction.bond_idx),
                            reaction.clone(),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .flatten()
            .collect();
    let molecules: Vec<MoleculeBlueprint> = reaction_registry
        .iter()
        .map(|(molecule, _)| molecule.clone())
        .collect();
    Some((forming_reactions, breaking_reactions, molecules))
}
