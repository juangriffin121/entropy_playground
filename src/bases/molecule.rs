use serde::Deserialize;
use std::{collections::HashMap, fs};

use super::recipient::{BlueprintAtomIndex, MoleculeAtomIndex, MoleculeId, WorldAtomId};

#[derive(Debug, Clone)]
pub struct Molecule {
    pub id: MoleculeId,
    pub atoms: Vec<WorldAtomId>,
    pub bonds_: Vec<(MoleculeAtomIndex, MoleculeAtomIndex)>, // (atom_1 pos in molecule, atom_2 pos in molecule)
    pub blueprint_id: usize,
}

#[derive(Debug, Deserialize)]
pub struct MoleculesJSON {
    pub blueprints: HashMap<String, MoleculeBlueprint>,
}
#[derive(Debug, Deserialize, Clone, PartialEq, Eq, Hash)]
pub struct MoleculeBlueprint {
    pub atoms: Vec<String>,
    pub bonds: Vec<(BlueprintAtomIndex, BlueprintAtomIndex)>, // Index of the first and second atom in the bond, index with respect to the
    // atoms vec in molecule blueprint
    pub canonical: bool,
}

impl MoleculesJSON {
    pub fn load(path: &str) -> Self {
        let json_content = fs::read_to_string(path).expect("Failed to read JSON file");
        serde_json::from_str(&json_content).expect("Failed to parse JSON")
    }
}

impl MoleculeBlueprint {
    pub fn canonicalize(
        &self,
    ) -> (
        MoleculeBlueprint,
        HashMap<BlueprintAtomIndex, BlueprintAtomIndex>,
    ) {
        let mut indices: Vec<usize> = (0..self.atoms.len()).collect();

        // Sort by atom type and number of bonds
        indices.sort_by_key(|&i| {
            (
                self.atoms[i].clone(),
                self.bonds
                    .iter()
                    .filter(|(a, b)| a.0 == i || b.0 == i)
                    .count(),
            )
        });

        // Reorder atoms and bonds
        let new_atoms: Vec<String> = indices.iter().map(|&i| self.atoms[i].clone()).collect();

        let mapping: HashMap<BlueprintAtomIndex, BlueprintAtomIndex> = indices
            .iter()
            .enumerate()
            .map(|(new_index, &old_index)| {
                (BlueprintAtomIndex(old_index), BlueprintAtomIndex(new_index))
            })
            .collect();

        let new_bonds = self
            .bonds
            .iter()
            .map(|(a, b)| {
                let new_a = indices.iter().position(|&x| x == a.0).unwrap();
                let new_b = indices.iter().position(|&x| x == b.0).unwrap();
                (
                    BlueprintAtomIndex(new_a.min(new_b)),
                    BlueprintAtomIndex(new_a.max(new_b)),
                )
            })
            .collect();
        // Create the canonicalized molecule
        let canonical_molecule = MoleculeBlueprint {
            atoms: new_atoms,
            bonds: new_bonds,
            canonical: true,
        };

        (canonical_molecule, mapping)
    }

    pub fn compare_canonical(&self, other: &MoleculeBlueprint) -> bool {
        if !(self.canonical && other.canonical) {
            panic!("tried to compare_canonical with non canonical molecules")
        }
        self.atoms == other.atoms && self.bonds == other.bonds
    }

    pub fn automatic_name(&self) -> String {
        let mut unique_count: HashMap<String, u32> = HashMap::new();
        println!("{:?}", self.atoms);
        for atom in &self.atoms {
            if unique_count.contains_key(atom) {
                *unique_count.get_mut(atom).unwrap() += 1
            } else {
                unique_count.insert(atom.to_string(), 1);
            }
        }

        let mut name = String::new();

        let mut sorted: Vec<(String, u32)> = unique_count.into_iter().collect();
        sorted.sort_by(|a, b| a.0.cmp(&b.0));
        for (atom, count) in sorted {
            name = format!("{}{}", name, atom);
            if count > 1 {
                name = format!("{}{}", name, count);
            }
        }

        name
    }
}
