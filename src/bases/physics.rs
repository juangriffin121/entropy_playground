use super::{
    atom::Atom,
    initialization::{Either, SingleOrPair},
    molecule::{Molecule, MoleculeBlueprint},
    reaction::{
        check_forming_reactions, FormingReactionBlueprint, FormingReactionKey, ParticleBlueprint,
    },
    recipient::{BlueprintAtomIndex, MoleculeId, Recipient, WorldAtomId},
};
use crate::bases::reaction::check_breaking_reactions;
use std::{
    collections::{HashMap, HashSet},
    fmt::Debug,
    iter::zip,
    usize,
};

pub enum Wall {
    Top,
    Bottom,
    Left,
    Right,
}

struct CollisionEngine;
impl CollisionEngine {
    pub fn check_collision(&self, a1: &Atom, a2: &Atom) -> bool {
        let dx = a1.position.0 - a2.position.0;
        let dy = a1.position.1 - a2.position.1;
        let r = a1.radius + a2.radius;
        dx.powi(2) + dy.powi(2) <= r.powi(2)
    }

    pub fn check_wall_collision(&self, atom: &Atom, shape: (f32, f32)) -> Option<Wall> {
        if atom.position.1 - atom.radius <= 0.0 {
            return Some(Wall::Top);
        } else if atom.position.1 + atom.radius >= shape.1 {
            return Some(Wall::Bottom);
        } else if atom.position.0 - atom.radius <= 0.0 {
            return Some(Wall::Left);
        } else if atom.position.0 + atom.radius >= shape.0 {
            return Some(Wall::Right);
        }
        None
    }

    pub fn process_collisions(&self, recipient: &mut Recipient) {
        let atoms = &mut recipient.contents;
        let dimensions = recipient.shape;
        let forming_reactions = &recipient.forming_reactions;
        let molecule_registry = &recipient.molecule_registry;
        let molecule_blueprints = &recipient.molecule_blueprints;
        let dt = recipient.dt;
        let reactions = match &mut recipient.grid {
            Some(grid) => {
                grid.populate(atoms);
                self.grid_based_collisions(
                    grid,
                    atoms,
                    dimensions,
                    forming_reactions,
                    molecule_registry,
                    molecule_blueprints,
                    dt,
                )
            }
            None => self.gridless_process_collisions(
                atoms,
                dimensions,
                forming_reactions,
                molecule_registry,
                molecule_blueprints,
                dt,
            ),
        };

        let molecule_registry = &mut recipient.molecule_registry;
        let bond_map = &mut recipient.bond_map;
        let next_molecule_id = &mut recipient.next_molecule_id;
        let blueprint_id_map = &recipient.blueprint_id_map;
        let bond_blueprints = &recipient.bond_blueprints;
        let mut count = 0;
        let mut molecules_seen: HashSet<MoleculeId> = HashSet::new();
        let mut atoms_seen: HashSet<WorldAtomId> = HashSet::new();

        for (atom_ids, reactant_idxs, reaction) in reactions {
            count += 1;
            println!("checking reaction {}", count);
            // Check if any reactant is already seen
            if any_in(&reactant_idxs, &mut molecules_seen) {
                println!("skipping reaction for seen reactant {:?}", reactant_idxs);
                continue;
            }
            if atoms_seen.contains(&atom_ids.0) || atoms_seen.contains(&atom_ids.1) {
                println!("skipping reaction for seen atom {:?}", atom_ids);
                continue;
            }
            atoms_seen.insert(atom_ids.0);
            atoms_seen.insert(atom_ids.1);
            reaction.apply(
                atom_ids,
                reactant_idxs,
                atoms,
                molecule_registry,
                bond_map,
                next_molecule_id,
                blueprint_id_map
                    .as_ref()
                    .expect("cant have reactions without molecule_blueprints"),
                bond_blueprints
                    .as_ref()
                    .expect("cant have reactions without bond_blueprints"),
            );
        }
    }

    pub fn grid_based_collisions<'a>(
        &self,
        grid: &Grid,
        atoms: &mut Vec<Atom>,
        dimensions: (f32, f32),
        forming_reactions: &'a Option<HashMap<FormingReactionKey, FormingReactionBlueprint>>,
        molecule_registry: &HashMap<MoleculeId, Molecule>,
        molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
        dt: f32,
    ) -> Vec<(
        (WorldAtomId, WorldAtomId),
        Either<
            MoleculeId,
            (
                Either<MoleculeId, WorldAtomId>,
                Either<MoleculeId, WorldAtomId>,
            ),
        >,
        &'a FormingReactionBlueprint,
    )> {
        let mut seen: HashSet<(WorldAtomId, WorldAtomId)> = HashSet::new();
        let mut reactions = Vec::new();
        for atom_idx in 0..atoms.len() {
            let cell = grid.get_cell(&atoms[atom_idx]);
            let neighbors = grid.neighbors(cell);
            for neighbor in neighbors {
                let candidate_atom_idxs = &grid.grid[neighbor.0][neighbor.1];
                for candidate_idx in candidate_atom_idxs {
                    if atom_idx == *candidate_idx {
                        continue;
                    }
                    let (i, j) = (
                        WorldAtomId(atom_idx.min(*candidate_idx)),
                        WorldAtomId(atom_idx.max(*candidate_idx)),
                    );

                    if seen.contains(&(i, j)) {
                        continue;
                    }
                    seen.insert((i, j));
                    let (left, right) = atoms.split_at_mut(j.0);
                    let a = &mut left[i.0];
                    let b = &mut right[0];
                    if self.check_collision(a, b) {
                        let reaction = self.resolve_particle_particle(
                            a,
                            b,
                            forming_reactions,
                            molecule_registry,
                            molecule_blueprints,
                            dt,
                        );
                        if let Some(reaction) = reaction {
                            reactions.push(reaction);
                        }
                    }
                }
            }

            if let Some(wall) = self.check_wall_collision(&atoms[atom_idx], dimensions) {
                self.resolve_particle_wall(&mut atoms[atom_idx], dimensions, wall);
            }
        }
        reactions
    }

    pub fn gridless_process_collisions<'a>(
        &self,
        atoms: &mut Vec<Atom>,
        dimensions: (f32, f32),
        forming_reactions: &'a Option<HashMap<FormingReactionKey, FormingReactionBlueprint>>,
        molecule_registry: &HashMap<MoleculeId, Molecule>,
        molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
        dt: f32,
    ) -> Vec<(
        (WorldAtomId, WorldAtomId),
        Either<
            MoleculeId,
            (
                Either<MoleculeId, WorldAtomId>,
                Either<MoleculeId, WorldAtomId>,
            ),
        >,
        &'a FormingReactionBlueprint,
    )> {
        let mut reactions = Vec::new();
        for i in 0..atoms.len() {
            for j in (i + 1)..atoms.len() {
                let (left, right) = atoms.split_at_mut(j);
                let a = &mut left[i];
                let b = &mut right[0];

                if self.check_collision(a, b) {
                    let reaction_info = self.resolve_particle_particle(
                        a,
                        b,
                        forming_reactions,
                        molecule_registry,
                        molecule_blueprints,
                        dt,
                    );
                    if let Some(reaction_info) = reaction_info {
                        reactions.push(reaction_info);
                    }
                }
            }

            if let Some(wall) = self.check_wall_collision(&atoms[i], dimensions) {
                self.resolve_particle_wall(&mut atoms[i], dimensions, wall);
            }
        }
        reactions
    }

    pub fn resolve_particle_particle<'a>(
        &self,
        a: &mut Atom,
        b: &mut Atom,
        forming_reactions: &'a Option<HashMap<FormingReactionKey, FormingReactionBlueprint>>,
        molecule_registry: &HashMap<MoleculeId, Molecule>,
        molecule_blueprints: &Option<Vec<MoleculeBlueprint>>,
        dt: f32,
    ) -> Option<(
        (WorldAtomId, WorldAtomId),
        Either<
            MoleculeId,
            (
                Either<MoleculeId, WorldAtomId>,
                Either<MoleculeId, WorldAtomId>,
            ),
        >,
        &'a FormingReactionBlueprint,
    )> {
        let dx = a.position.0 - b.position.0;
        let dy = a.position.1 - b.position.1;
        let distance = (dx * dx + dy * dy).sqrt();

        if distance < 1e-6 {
            println!("happened {:?}, {:?}", a.id, b.id);
            return None;
        }

        let nx = dx / distance;
        let ny = dy / distance;

        let dvx = a.velocity.0 - b.velocity.0;
        let dvy = a.velocity.1 - b.velocity.1;

        let dot_product = dvx * nx + dvy * ny;

        if dt * dot_product >= 1e-6 {
            return None;
        }

        let mass1 = a.mass;
        let mass2 = b.mass;
        let mass_sum = mass1 + mass2;

        let impulse = 2.0 * dot_product / mass_sum;

        let reaction_data = check_forming_reactions(
            a,
            b,
            forming_reactions,
            molecule_registry,
            molecule_blueprints,
        );

        a.velocity.0 -= impulse * mass2 * nx;
        a.velocity.1 -= impulse * mass2 * ny;
        b.velocity.0 += impulse * mass1 * nx;
        b.velocity.1 += impulse * mass1 * ny;

        if let Some((idxs, reaction)) = reaction_data {
            Some(((a.id, b.id), idxs, reaction))
        } else {
            None
        }
    }

    pub fn resolve_particle_wall(&self, atom: &mut Atom, shape: (f32, f32), wall: Wall) {
        match wall {
            Wall::Top => {
                atom.position.1 = atom.radius - (atom.position.1 - atom.radius);
                atom.velocity.1 = -atom.velocity.1;
            }
            Wall::Bottom => {
                atom.position.1 = shape.1 - atom.radius - (atom.position.1 + atom.radius - shape.1);
                atom.velocity.1 = -atom.velocity.1;
            }
            Wall::Left => {
                atom.position.0 = atom.radius - (atom.position.0 - atom.radius);
                atom.velocity.0 = -atom.velocity.0;
            }
            Wall::Right => {
                atom.position.0 = shape.0 - atom.radius - (atom.position.0 + atom.radius - shape.0);
                atom.velocity.0 = -atom.velocity.0;
            }
        }
    }
}

// Possible optimization
// Using Vec<Vec<Vec<usize>>> may lead to high memory usage if the grid is large but sparsely populated.
// Consider using a flat Vec<Vec<usize>> indexed by (x * grid_width + y) for potentially better cache locality.
pub struct Grid {
    pub grid: Vec<Vec<Vec<usize>>>, //position in grid -> Vec of positions in contents Vec
    pub cells: (i32, i32),
    pub cell_size: (f32, f32),
}

#[derive(Clone, serde::Serialize, serde::Deserialize, Debug)]
pub struct GridJSON {
    pub cells: (i32, i32),
}

impl Grid {
    pub fn new(cells: (i32, i32), cell_size: (f32, f32)) -> Self {
        let grid = vec![vec![Vec::new(); cells.1 as usize]; cells.0 as usize];
        Grid {
            grid,
            cells,
            cell_size,
        }
    }

    pub fn populate(&mut self, atoms: &Vec<Atom>) {
        self.clear();
        for (i, atom) in atoms.iter().enumerate() {
            self.place(atom, i as usize);
        }
    }

    pub fn from_grid_json(data: GridJSON, recipient_shape: (f32, f32)) -> Grid {
        let cells = data.cells;
        let cell_size = (
            recipient_shape.0 / cells.0 as f32,
            recipient_shape.1 / cells.1 as f32,
        );
        Grid::new(cells, cell_size)
    }

    pub fn get_cell(&self, atom: &Atom) -> (usize, usize) {
        let x = ((atom.position.0 / self.cell_size.0).floor() as isize).max(0) as usize;
        let y = ((atom.position.1 / self.cell_size.1).floor() as isize).max(0) as usize;
        (
            x.min(self.cells.0 as usize - 1),
            y.min(self.cells.1 as usize - 1),
        )
    }

    pub fn place(&mut self, atom: &Atom, index: usize) {
        let cell = self.get_cell(atom);
        self.grid[cell.0][cell.1].push(index);
    }

    pub fn clear(&mut self) {
        for row in &mut self.grid {
            for cell in row {
                cell.clear();
            }
        }
    }
    pub fn neighbors(&self, cell: (usize, usize)) -> Vec<(usize, usize)> {
        let mut neighbors = Vec::new();

        let (x, y) = cell;
        let (max_x, max_y) = (self.cells.0 as usize, self.cells.1 as usize);

        for dx in -1..=1 {
            for dy in -1..=1 {
                let nx = x as isize + dx;
                let ny = y as isize + dy;

                if nx >= 0 && nx < max_x as isize && ny >= 0 && ny < max_y as isize {
                    neighbors.push((nx as usize, ny as usize));
                }
            }
        }

        neighbors
    }
}

pub trait Force: Send + Sync {
    fn calculate(&self, atom: &Atom, recipient: &Recipient) -> (f32, f32);
}

pub struct Gravity {
    pub g: f64,
}

impl Force for Gravity {
    fn calculate(&self, atom: &Atom, recipient: &Recipient) -> (f32, f32) {
        return (0.0, self.g as f32);
    }
}

impl Atom {
    fn update(&mut self, dt: f32, force: (f32, f32)) {
        self.velocity.0 += force.0 * dt;
        self.velocity.1 += force.1 * dt;
        self.position.0 += self.velocity.0 * dt;
        self.position.1 += self.velocity.1 * dt;
    }
}

impl Recipient {
    pub fn update(&mut self, dt: f32) {
        self.update_bonds();

        let forces: Vec<(f32, f32)> = self
            .contents
            .iter()
            .map(|atom| self.total_force(atom))
            .collect();
        for (atom, force) in zip(&mut self.contents, forces) {
            atom.update(dt, force);
        }

        CollisionEngine.process_collisions(self);
        println!("molecules: {}", self.molecule_registry.len());
    }

    pub fn total_force(&self, atom: &Atom) -> (f32, f32) {
        self.force_registry
            .iter()
            .map(|force| force.calculate(atom, self))
            .fold((0.0, 0.0), |(acc_x, acc_y), (f_x, f_y)| {
                (acc_x + f_x, acc_y + f_y)
            })
    }

    fn update_bonds(&mut self) {
        // First collect all updates without mutating anything
        let bond_updates: Vec<((WorldAtomId, WorldAtomId), (f32, f32), bool)> = self
            .bond_map
            .iter()
            .map(|(&(atom1_id, atom2_id), bond)| {
                let atom1 = &self.contents[atom1_id.0];
                let atom2 = &self.contents[atom2_id.0];
                let dx = atom2.position.0 - atom1.position.0;
                let dy = atom2.position.1 - atom1.position.1;
                let distance = (dx * dx + dy * dy).sqrt();
                let force_magnitude =
                    bond.k * (distance - (bond.equilibrium_distance + atom1.radius + atom2.radius));
                let fx = force_magnitude * dx / distance;
                let fy = force_magnitude * dy / distance;
                (
                    (atom1_id, atom2_id),
                    (fx, fy),
                    distance > bond.breaking_distance + atom1.radius + atom2.radius,
                )
            })
            .collect();

        // Apply updates and collect bonds to break
        let mut broken_bonds = Vec::new();
        for ((atom1_id, atom2_id), force, should_break) in bond_updates {
            if let Some(bond) = self.bond_map.get_mut(&(atom1_id, atom2_id)) {
                bond.force = force;
            }
            if should_break {
                broken_bonds.push((atom1_id, atom2_id));
            }
        }

        // Remove broken bonds REACTIONNNNNNNNN
        for (atom1_id, atom2_id) in broken_bonds {
            let atom1 = &self.contents[atom1_id.0];
            let atom2 = &self.contents[atom2_id.0];
            let reaction = check_breaking_reactions(
                atom1,
                atom2,
                &self.breaking_reactions,
                &self.molecule_registry,
                &self.molecule_blueprints,
            );
            if let Some((molecule_id, reaction)) = reaction {
                reaction.apply(
                    atom1_id,
                    atom2_id,
                    &mut self.contents,
                    molecule_id,
                    &mut self.molecule_registry,
                    &mut self.bond_map,
                    &mut self.next_molecule_id,
                    self.molecule_blueprints.as_ref().unwrap(),
                    self.blueprint_id_map.as_ref().unwrap(),
                );
            }
        }
    }
}

fn any_in(
    reactant_idxs: &Either<
        MoleculeId,
        (
            Either<MoleculeId, WorldAtomId>,
            Either<MoleculeId, WorldAtomId>,
        ),
    >,
    seen: &mut HashSet<MoleculeId>,
) -> bool {
    match reactant_idxs {
        Either::Left(mol_id) => {
            let answer = seen.contains(&mol_id);
            if !answer {
                seen.insert(mol_id.clone());
            }
            answer
        }
        Either::Right((id1, id2)) => {
            let answer1 = match id1 {
                Either::Left(mol_id) => {
                    let answer = seen.contains(&mol_id);
                    if !answer {
                        seen.insert(mol_id.clone());
                    }
                    answer
                }
                Either::Right(_) => false,
            };
            let answer2 = match id2 {
                Either::Left(mol_id) => {
                    let answer = seen.contains(&mol_id);
                    if !answer {
                        seen.insert(mol_id.clone());
                    }
                    answer
                }
                Either::Right(_) => false,
            };
            answer1 || answer2
        }
    }
}
